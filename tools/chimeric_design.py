#!/usr/bin/env python3
"""
Chimeric Protein Design Pipeline for PH20M3
============================================

Steps:
  1. Read PH20M3 sequence and MSA A3M file
  2. Define loops to be replaced (e.g., "270-284;68-74")
  3. Find homologous loops in MSA sequences (Homo sapiens, overall identity > 0.3)
  4. Build chimeric sequences by swapping loops
  5. Create Chai-1 SLURM prediction tasks
  6. Convert model_idx_0.cif → .pdb for each chimeric sequence
  7. Score PDBs with ProteinMPNN
  8. argparse for input/output control

Usage:
  python chimeric_design.py \\
    --fasta PH20M3.fasta \\
    --msa PH20M3_uniprot.a3m \\
    --loops "270-284;68-74" \\
    --output ./chimeric_results \\
    --min-identity 0.3 \\
    --max-chimeras 50 \\
    --submit-slurm \\
    --slurm-partition 4090 \\
    --ncpus 4
"""

import argparse
import os
import sys
import re
import subprocess as sp
import shutil
import textwrap
import csv
from collections import OrderedDict
from itertools import product as cartprod


# ---------------------------------------------------------------------------
# Paths (adjust if needed)
# ---------------------------------------------------------------------------
CHAI1_DIR = "/data_test/share/pub_tools/chai-lab"
CHAI1_RUN = os.path.join(CHAI1_DIR, "run.sh")
SLURM_SUBMIT = "/data_test/home/lzzheng/bin/submit_slurm_gpu.sh"
PROTEINMPNN_DIR = "/data_test/home/lzzheng/apps/ProteinMPNN"
PROTEINMPNN_WRAPPER = os.path.join(PROTEINMPNN_DIR, "run_ProteinMPNN.py")
CIF2PDB_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_cif2pdb.py")


# ---------------------------------------------------------------------------
# 1.  FASTA / A3M readers
# ---------------------------------------------------------------------------
def read_fasta(fpath):
    """Return (header, sequence) for a single-record FASTA."""
    header = None
    lines = []
    with open(fpath) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    break
                header = line
            else:
                lines.append(line)
    if header is None:
        raise ValueError(f"No FASTA header found in {fpath}")
    return header, "".join(lines)


def parse_a3m(fpath):
    """Return list of (header, sequence) from an A3M file."""
    entries = []
    curr_hdr = None
    curr_seq = []
    with open(fpath) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if curr_hdr is not None:
                    entries.append((curr_hdr, "".join(curr_seq)))
                curr_hdr = line
                curr_seq = []
            else:
                curr_seq.append(line)
    if curr_hdr is not None:
        entries.append((curr_hdr, "".join(curr_seq)))
    return entries


# ---------------------------------------------------------------------------
# 2.  Loop definition helpers
# ---------------------------------------------------------------------------
def parse_loops(loop_str):
    """Parse e.g. '270-284;68-74' → [(270,284), (68,74)] (1-based inclusive)."""
    pairs = []
    for part in loop_str.replace(";", ",").split(","):
        part = part.strip()
        if not part:
            continue
        m = re.match(r"^(\d+)\s*-\s*(\d+)$", part)
        if not m:
            raise ValueError(f"Cannot parse loop spec: '{part}'")
        start, end = int(m.group(1)), int(m.group(2))
        if start > end:
            start, end = end, start
        pairs.append((start, end))
    return pairs


# ---------------------------------------------------------------------------
# 3.  Sequence identity & human filter
# ---------------------------------------------------------------------------
def seq_identity(seq_a, seq_b):
    """Fraction of identical residues over aligned positions (ignoring gaps in seq_b)."""
    matches = total = 0
    for a, b in zip(seq_a, seq_b):
        if b == "-":
            continue
        total += 1
        if a == b:
            matches += 1
    return matches / total if total else 0.0


def extract_loop_region(a3m_seq, query_seq, q_start, q_end):
    """
    Given an A3M-aligned sequence and the query (PH20M3), extract the
    sub-sequence corresponding to query positions q_start..q_end (1-based).

    Walks through the alignment counting non-gap positions in the query.
    Returns the donor sequence with gaps stripped.
    """
    q_idx = 0  # number of non-gap positions seen in query
    donor_chars = []
    for q_char, d_char in zip(query_seq, a3m_seq):
        if q_char != "-":
            q_idx += 1
        if q_start <= q_idx <= q_end:
            if d_char != "-":
                donor_chars.append(d_char.upper())
    return "".join(donor_chars)


def find_human_homologs(msa_entries, query_seq, min_identity=0.3):
    """
    Filter MSA entries to Homo sapiens sequences with overall identity > min_identity.

    Returns list of (header, aligned_seq, identity, short_label).
    """
    results = []
    for hdr, aln_seq in msa_entries:
        if "Homo sapiens" not in hdr and "_HUMAN" not in hdr:
            continue
        # Use the full-length aligned sequence for identity
        ident = seq_identity(query_seq, aln_seq)
        if ident < min_identity:
            continue
        # Build a short label from the header
        label = hdr.split()[0].lstrip(">").split("|")[-1] if "|" in hdr else hdr.split()[0].lstrip(">")
        results.append((hdr, aln_seq, ident, label))
    return results


def sanitize_name(name):
    """Replace characters unsafe for filenames."""
    return re.sub(r'[/\\:<>"|?* ]', "_", name)


# ---------------------------------------------------------------------------
# 4.  Chimeric sequence builder
# ---------------------------------------------------------------------------
def build_chimeras(query_seq, homologs, loops, max_chimeras=200):
    """
    Generate ALL combinatorial chimeras by replacing each loop region with
    every unique donor from human homologs.

    Returns:
        chimeras: list of (name, sequence)
        chimera_loops: dict mapping name -> list of (q_start, q_end, donor_label, donor_seq)
    """
    query_len = len(query_seq)
    chimera_loops = OrderedDict()

    # Collect donor options per loop (deduplicated by sequence)
    loop_options = OrderedDict()
    for q_start, q_end in loops:
        seen_seq = OrderedDict()
        for hdr, aln_seq, ident, label in homologs:
            donor_loop = extract_loop_region(aln_seq, query_seq, q_start, q_end)
            if len(donor_loop) < 2:
                continue
            if donor_loop not in seen_seq:
                seen_seq[donor_loop] = (label, ident)
        # Always include wild-type as an option (the original PH20M3 loop)
        wt_loop = query_seq[q_start - 1 : q_end]
        if wt_loop not in seen_seq:
            seen_seq[wt_loop] = ("PH20M3_wt", 1.0)
        # Move wild-type to front so combo[0] is always a donor (skip wt-only)
        wt_key = None
        for k in seen_seq:
            if k == wt_loop:
                wt_key = k
                break
        if wt_key:
            seen_seq.move_to_end(wt_key, last=False)
        loop_options[(q_start, q_end)] = [
            (label, dl, ident) for dl, (label, ident) in seen_seq.items()
        ]

    # Report available donor loops
    total_combos = 1
    for (q_start, q_end), opts in loop_options.items():
        n_donors = len(opts) - 1  # exclude wild-type
        print(f"  Loop {q_start}-{q_end} ({query_seq[q_start-1:q_end]}): "
              f"{n_donors} non-WT donors, {len(opts)} total options")
        for label, donor_loop, ident in opts:
            marker = " [WT]" if donor_loop == query_seq[q_start-1:q_end] else ""
            print(f"    {label}{marker} (id={ident:.3f}) {donor_loop}")
        total_combos *= len(opts)

    print(f"  Total possible combinations (incl. WT): {total_combos}")

    # Build ALL combinations via Cartesian product
    all_opts = list(loop_options.values())
    loop_keys = list(loop_options.keys())
    combo_idx = 0
    for combo in cartprod(*all_opts):
        # combo is a tuple of (label, donor_loop, ident) for each loop
        # Skip the all-wild-type combination (it's just PH20M3)
        if all(dl == query_seq[s - 1 : e] for (s, e), (_, dl, _) in zip(loop_keys, combo)):
            continue

        # Build the chimeric sequence
        seq_list = list(query_seq)
        loop_info = []
        name_parts = []
        for (q_start, q_end), (label, donor_loop, ident) in zip(loops, combo):
            if donor_loop != query_seq[q_start - 1 : q_end]:
                seq_list[q_start - 1 : q_end] = list(donor_loop)
                name_parts.append(f"L{q_start}-{q_end}_{label}")
            loop_info.append((q_start, q_end, label, donor_loop))
        new_seq = "".join(seq_list)

        name = "chimera_" + "_".join(name_parts) if name_parts else "chimera_all_wt"
        chimera_loops[name] = loop_info
        combo_idx += 1

    if combo_idx > max_chimeras:
        print(f"  Limiting to {max_chimeras} chimeras (from {combo_idx})")
        # Keep a diverse subset: take from middle (skip first few all-WT-like)
        keys = list(chimera_loops.keys())
        step = len(keys) // max_chimeras
        selected_keys = keys[1::step][:max_chimeras]  # skip first (WT-only)
        chimera_loops = OrderedDict((k, chimera_loops[k]) for k in selected_keys)

    chimeras = [(name, build_chimera_seq(query_seq, loops, chimera_loops[name]))
                for name in chimera_loops]
    return chimeras, chimera_loops


def build_chimera_seq(query_seq, loops, loop_info):
    """Apply loop replacements to query_seq.  loop_info is list of
    (q_start, q_end, label, donor_seq) — only entries where donor != WT
    have been included."""
    seq_list = list(query_seq)
    for q_start, q_end, label, donor_seq in loop_info:
        if donor_seq != query_seq[q_start - 1 : q_end]:
            seq_list[q_start - 1 : q_end] = list(donor_seq)
    return "".join(seq_list)


# ---------------------------------------------------------------------------
# 5.  Write FASTA for Chai-1
# ---------------------------------------------------------------------------
def write_fasta(seq, outpath, header):
    """Write a single-record FASTA with 'protein|' prefix for Chai-1 compatibility."""
    with open(outpath, "w") as fh:
        fh.write(f">protein|{header}\n")
        # wrap at 60 chars
        for i in range(0, len(seq), 60):
            fh.write(seq[i : i + 60] + "\n")


# ---------------------------------------------------------------------------
# 6.  CIF → PDB conversion (via BioPython or simple helper)
# ---------------------------------------------------------------------------
def make_cif2pdb_script():
    """Write a small helper script that converts a CIF file to PDB."""
    code = '''#!/usr/bin/env python3
"""Convert Chai-1 CIF output to PDB using BioPython."""
import sys
try:
    from Bio.PDB import MMCIFParser, PDBIO
except ImportError:
    print("BioPython not available; skipping CIF->PDB", file=sys.stderr)
    sys.exit(1)

cif_path = sys.argv[1]
pdb_path = sys.argv[2]
parser = MMCIFParser()
structure = parser.get_structure("model", cif_path)
io = PDBIO()
io.set_structure(structure)
io.save(pdb_path)
print(f"Written {pdb_path}")
'''
    with open(CIF2PDB_SCRIPT, "w") as fh:
        fh.write(code)
    os.chmod(CIF2PDB_SCRIPT, 0o755)
    return CIF2PDB_SCRIPT


def convert_cif_to_pdb(cif_path, pdb_path):
    """Run the external cif→pdb script."""
    if not os.path.exists(CIF2PDB_SCRIPT):
        make_cif2pdb_script()
    result = sp.run([sys.executable, CIF2PDB_SCRIPT, cif_path, pdb_path],
                    stdout=sp.PIPE, stderr=sp.PIPE)
    if result.returncode != 0:
        print(f"  WARNING: cif→pdb failed: {result.stderr.decode().strip()}", file=sys.stderr)
        return False
    return True


# ---------------------------------------------------------------------------
# 7.  ProteinMPNN scoring
# ---------------------------------------------------------------------------
def score_pdb_with_proteinmpnn(pdb_path, output_dir, seq_name, chimera_seq=None):
    """Run ProteinMPNN scoring on a single PDB; returns path to output FASTA.

    Uses the existing run_score_only.sh script which works with the sfct environment.
    """
    mpnn_out = os.path.join(output_dir, f"mpnn_{seq_name}")
    os.makedirs(mpnn_out, exist_ok=True)

    # Use the run_score_only.sh script
    run_score_script = "/data_test/home/lzzheng/apps/ProteinMPNN/run_score_only.sh"
    if not os.path.exists(run_score_script):
        print(f"  ERROR: run_score_only.sh not found at {run_score_script}", file=sys.stderr)
        return None

    # Create a FASTA file with the chimeric sequence
    fasta_path = os.path.join(mpnn_out, f"{seq_name}.fasta")
    if chimera_seq:
        # Use provided sequence
        with open(fasta_path, "w") as fh:
            fh.write(f">protein|{seq_name}\n{chimera_seq}\n")
    else:
        # Try to extract sequence from PDB (fallback)
        # For now, just create a dummy FASTA - we need the actual sequence
        print(f"  WARNING: No sequence provided for {seq_name}, using dummy", file=sys.stderr)
        with open(fasta_path, "w") as fh:
            fh.write(f">protein|{seq_name}\nM\n")  # Dummy

    # Add chain ID (default to "A")
    cmd = [run_score_script, pdb_path, fasta_path, mpnn_out, "A"]
    print(f"  Running ProteinMPNN: {run_score_script} {pdb_path} {fasta_path} {mpnn_out} A")

    # Suppress warnings
    env = os.environ.copy()
    env["PYTHONWARNINGS"] = "ignore"

    result = sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE, env=env)
    if result.returncode != 0:
        stderr_text = result.stderr.decode()
        if "Warning: importing 'simtk.openmm' is deprecated" in stderr_text:
            # Just a warning, check if output was generated
            pass
        else:
            print(f"  ProteinMPNN failed (exit {result.returncode}): {stderr_text[:200]}", file=sys.stderr)

    # Find output NPZ files (scores)
    npz_files = [f for f in os.listdir(mpnn_out) if f.endswith(".npz")]
    # Also check score_only subdirectory
    score_only_dir = os.path.join(mpnn_out, "score_only")
    if os.path.isdir(score_only_dir):
        npz_files.extend([os.path.join("score_only", f) for f in os.listdir(score_only_dir) if f.endswith(".npz")])

    if npz_files:
        # Return first NPZ file path
        return os.path.join(mpnn_out, npz_files[0])

    # Fallback: look for FASTA output (older style)
    seqs_dir = os.path.join(mpnn_out, "seqs")
    if os.path.isdir(seqs_dir):
        fas = [f for f in os.listdir(seqs_dir) if f.endswith(".fa")]
        if fas:
            return os.path.join(seqs_dir, fas[0])
    fas = [f for f in os.listdir(mpnn_out) if f.endswith(".fa")]
    if fas:
        return os.path.join(mpnn_out, fas[0])

    return None


def extract_mpnn_scores(fasta_path):
    """Parse ProteinMPNN FASTA to extract scores from headers."""
    scores = {}
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                m = re.search(r"score=([-\d.]+)", line)
                if m:
                    scores[line.strip()] = float(m.group(1))
    return scores


def extract_mpnn_scores_from_npz(npz_path):
    """Extract best score from ProteinMPNN NPZ file."""
    try:
        import numpy as np
    except ImportError:
        print(f"  numpy not available, cannot read {npz_path}", file=sys.stderr)
        return None

    try:
        data = np.load(npz_path)
        # Look for score arrays
        if "score" in data:
            scores = data["score"]
            if scores.size > 0:
                best = float(np.min(scores))
                data.close()
                return best
        elif "global_score" in data:
            scores = data["global_score"]
            if scores.size > 0:
                best = float(np.min(scores))
                data.close()
                return best
        data.close()
    except Exception as e:
        print(f"  Error reading NPZ {npz_path}: {e}", file=sys.stderr)
    return None


# ---------------------------------------------------------------------------
# 8.  Chai-1 submission
# ---------------------------------------------------------------------------
def submit_chai1(fasta_path, output_dir, msa_path=None, submit_slurm=False,
                 slurm_partition="4090", ncpus=4):
    """Run or submit Chai-1 prediction."""
    if submit_slurm:
        # Build the command string for SLURM
        chai_cmd = f"{CHAI1_RUN} {fasta_path} {output_dir}"
        if msa_path:
            chai_cmd += f" {msa_path}"
        slurm_cmd = f"{SLURM_SUBMIT} \"{chai_cmd}\" {ncpus} {slurm_partition}"
        print(f"  Submitting SLURM job: {slurm_cmd}")
        result = sp.run(slurm_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        out = (result.stdout or b"").decode() + (result.stderr or b"").decode()
        print(f"  SLURM output: {out.strip()}")
        return out
    else:
        # Run directly
        cmd = [CHAI1_RUN, fasta_path, output_dir]
        if msa_path:
            cmd.append(msa_path)
        print(f"  Running: {' '.join(cmd)}")
        result = sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        if result.returncode != 0:
            print(f"  WARNING: Chai-1 failed: {result.stderr.decode()[:200]}", file=sys.stderr)
        return (result.stdout or b"").decode() or (result.stderr or b"").decode()


# ---------------------------------------------------------------------------
# 9.  Score summary
# ---------------------------------------------------------------------------
# Chai-1 score extraction from NPZ files
# ---------------------------------------------------------------------------
def load_chai_scores(npz_path):
    """Load Chai-1 NPZ scores.  Returns dict with plddt, ptm, iptm, etc."""
    try:
        import numpy as np
    except ImportError:
        print("  numpy not available; cannot read NPZ scores", file=sys.stderr)
        return {}
    try:
        data = np.load(npz_path)
        scores = {}
        for key in data.keys():
            val = data[key]
            if val.ndim == 0:
                scores[key] = float(val)
            elif val.ndim == 1 and val.size == 1:
                scores[key] = float(val[0])
            elif key == "plddt" and val.ndim == 1:
                scores["plddt_mean"] = float(val.mean())
                scores["plddt"] = val.tolist()
            else:
                scores[key] = val.tolist()
        data.close()
        return scores
    except Exception as e:
        print(f"  Error loading NPZ {npz_path}: {e}", file=sys.stderr)
        return {}


def extract_chai_confidence(cif_path):
    """Extract average pLDDT from a CIF file (B-factor column)."""
    if not os.path.exists(cif_path):
        return None

    # Try BioPython first
    try:
        from Bio.PDB import MMCIFParser
        parser = MMCIFParser()
        structure = parser.get_structure("model", cif_path)
        b_factors = []
        for atom in structure.get_atoms():
            if atom.get_name() == "CA":
                b_factors.append(atom.get_bfactor())
        return sum(b_factors) / len(b_factors) if b_factors else None
    except ImportError:
        pass  # Fall back to manual parsing
    except Exception:
        pass  # Fall back to manual parsing

    # Manual CIF parsing for B-factors (CA atoms only)
    try:
        with open(cif_path) as fh:
            lines = fh.readlines()

        # Find the atom_site loop and column indices
        in_atom_site = False
        b_factor_idx = -1
        atom_name_idx = -1
        col_idx = 0
        b_factors = []

        for line in lines:
            line = line.strip()
            if line.startswith("loop_"):
                in_atom_site = False
                b_factor_idx = -1
                atom_name_idx = -1
                col_idx = 0
            elif line.startswith("_atom_site."):
                in_atom_site = True
                if "B_iso_or_equiv" in line:
                    b_factor_idx = col_idx
                elif "label_atom_id" in line:  # Use label_atom_id not auth_atom_id
                    atom_name_idx = col_idx
                col_idx += 1
            elif in_atom_site and line and not line.startswith("#") and not line.startswith("_"):
                parts = line.split()
                if len(parts) > max(b_factor_idx, atom_name_idx) >= 0:
                    atom_name = parts[atom_name_idx]
                    if atom_name == "CA" and b_factor_idx >= 0:
                        try:
                            b_factors.append(float(parts[b_factor_idx]))
                        except ValueError:
                            pass

        return sum(b_factors) / len(b_factors) if b_factors else None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# CSV report
# ---------------------------------------------------------------------------
def write_full_report(report_path, chimeras, chimera_loops, query_seq, loops,
                      chai_scores, mpnn_scores):
    """Write a comprehensive CSV with loop selections, sequence, scores."""
    fieldnames = [
        "name", "sequence", "length",
    ]
    # Add loop donor columns
    for q_start, q_end in loops:
        fieldnames.append(f"loop_{q_start}-{q_end}_donor")
        fieldnames.append(f"loop_{q_start}-{q_end}_seq")
    fieldnames += [
        "chai1_plddt", "chai1_ptm", "chai1_iptm",
        "mpnn_best_score", "rank_combined",
    ]

    rows = []
    for name, seq in chimeras:
        safe = sanitize_name(name)
        row = {"name": name, "sequence": seq, "length": len(seq)}

        # Loop info
        loop_info_dict = {(s, e): (lbl, ds) for s, e, lbl, ds in chimera_loops.get(name, [])}
        for q_start, q_end in loops:
            if (q_start, q_end) in loop_info_dict:
                lbl, ds = loop_info_dict[(q_start, q_end)]
                row[f"loop_{q_start}-{q_end}_donor"] = lbl
                row[f"loop_{q_start}-{q_end}_seq"] = ds
            else:
                row[f"loop_{q_start}-{q_end}_donor"] = "PH20M3_wt"
                row[f"loop_{q_start}-{q_end}_seq"] = query_seq[q_start-1:q_end]

        # Chai-1 scores
        cs = chai_scores.get(safe, {})
        row["chai1_plddt"] = f"{cs.get('plddt_mean', ''):.2f}" if "plddt_mean" in cs else ""
        row["chai1_ptm"] = f"{cs.get('ptm', ''):.4f}" if "ptm" in cs else ""
        row["chai1_iptm"] = f"{cs.get('iptm', ''):.4f}" if "iptm" in cs else ""

        # MPNN score
        row["mpnn_best_score"] = f"{mpnn_scores.get(safe, ''):.4f}" if safe in mpnn_scores else ""

        rows.append(row)

    # Rank by combined score (pLDDT + pTM, favoring high pLDDT)
    def combined_rank(row):
        plddt = float(row["chai1_plddt"]) if row["chai1_plddt"] else 0
        ptm = float(row["chai1_ptm"]) if row["chai1_ptm"] else 0
        mpnn = float(row["mpnn_best_score"]) if row["mpnn_best_score"] else 0
        # Higher pLDDT/pTM is better, lower MPNN is better
        return -(plddt * 0.7 + ptm * 30) + mpnn * 0.3

    rows.sort(key=combined_rank)
    for i, row in enumerate(rows, 1):
        row["rank_combined"] = i

    with open(report_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Full report written to {report_path}")
    print(f"  Rows: {len(rows)}")
    if rows:
        print(f"  Top chimera: {rows[0]['name']} "
              f"(pLDDT={rows[0]['chai1_plddt']}, pTM={rows[0]['chai1_ptm']}, "
              f"MPNN={rows[0]['mpnn_best_score']})")
    return rows


def write_score_summary(scores, outpath):
    """Write a CSV summary of MPNN scores."""
    with open(outpath, "w") as fh:
        fh.write("seq_name,mpnn_score\n")
        for name, score in sorted(scores.items(), key=lambda x: x[1]):
            fh.write(f"{name},{score:.4f}\n")
    print(f"  Score summary written to {outpath}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Chimeric PH20M3 design pipeline: find homologous loops in "
                    "human MSA, build chimeras, predict with Chai-1, score with MPNN.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples:
              # Dry-run: generate chimeras and FASTA files only
              %(prog)s --fasta PH20M3.fasta --msa PH20M3_uniprot.a3m \\
                       --loops "270-284;68-74" --output ./chimeric_test

              # Full pipeline with SLURM submission
              %(prog)s --fasta PH20M3.fasta --msa PH20M3_uniprot.a3m \\
                       --loops "270-284;68-74" --output ./chimeric_full \\
                       --submit-slurm --slurm-partition 4090 --ncpus 4

              # Score existing PDBs (skip Chai-1)
              %(prog)s --fasta PH20M3.fasta --msa PH20M3_uniprot.a3m \\
                       --loops "270-284;68-74" --output ./chimeric_full \\
                       --score-only
        """),
    )
    ap.add_argument("--fasta", required=True, help="PH20M3 FASTA file")
    ap.add_argument("--msa", required=True, help="MSA A3M file (UniProt search result)")
    ap.add_argument("--loops", required=True,
                    help="Loop regions to swap, e.g. '270-284;68-74'")
    ap.add_argument("--output", "-o", default="./chimeric_results",
                    help="Output directory (default: ./chimeric_results)")
    ap.add_argument("--min-identity", type=float, default=0.3,
                    help="Minimum overall sequence identity to human homolog (default: 0.3)")
    ap.add_argument("--max-chimeras", type=int, default=200,
                    help="Maximum number of chimeric sequences to generate (default: 200)")
    ap.add_argument("--submit-slurm", action="store_true",
                    help="Submit Chai-1 jobs via SLURM instead of running directly")
    ap.add_argument("--slurm-partition", default="4090",
                    help="SLURM partition (default: 4090)")
    ap.add_argument("--ncpus", type=int, default=4,
                    help="Number of CPUs for SLURM job (default: 4)")
    ap.add_argument("--score-only", action="store_true",
                    help="Skip Chai-1 prediction; only run MPNN scoring on existing PDBs")
    ap.add_argument("--report-only", action="store_true",
                    help="Collate scores from existing predictions into CSV report (no new jobs)")
    ap.add_argument("--skip-chai1", action="store_true",
                    help="Skip Chai-1 prediction entirely (only generate chimeras)")
    ap.add_argument("--skip-mpnn", action="store_true",
                    help="Skip ProteinMPNN scoring")
    ap.add_argument("--msa-for-chai1",
                    help="Optional A3M MSA to pass to Chai-1 predictions")

    args = ap.parse_args()

    # -----------------------------------------------------------------------
    # Setup output directory structure
    # -----------------------------------------------------------------------
    root = os.path.abspath(args.output)
    fasta_dir = os.path.join(root, "fastas")
    chai_dir = os.path.join(root, "chai1_results")
    pdb_dir = os.path.join(root, "pdbs")
    mpnn_dir = os.path.join(root, "mpnn_scores")
    os.makedirs(fasta_dir, exist_ok=True)
    os.makedirs(chai_dir, exist_ok=True)
    os.makedirs(pdb_dir, exist_ok=True)
    os.makedirs(mpnn_dir, exist_ok=True)

    # Fast-path: report-only mode — collate from existing results
    if args.report_only:
        print("\n" + "=" * 60)
        print("REPORT-ONLY MODE: collating scores from existing predictions")
        print("=" * 60)

        # Read existing chimeras index
        index_csv = os.path.join(root, "chimeras_index.csv")
        if not os.path.exists(index_csv):
            print(f"  ERROR: chimera index not found at {index_csv}", file=sys.stderr)
            sys.exit(1)
        chimeras = []
        chimera_loops = {}
        loops = parse_loops(args.loops)
        _, query_seq = read_fasta(args.fasta)
        with open(index_csv) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                name = row["name"]
                loop_info = []
                for q_start, q_end in loops:
                    donor = row.get(f"loop_{q_start}-{q_end}_donor", "PH20M3_wt")
                    dseq = row.get(f"loop_{q_start}-{q_end}_seq", query_seq[q_start-1:q_end])
                    loop_info.append((q_start, q_end, donor, dseq))
                chimera_loops[name] = loop_info
                seq = row.get("sequence", "")
                chimeras.append((name, seq))

        print(f"  Loaded {len(chimeras)} chimeras from index")

        # Extract Chai-1 scores and MPNN scores
        chai_scores = {}
        for name, _ in chimeras:
            safe_name = sanitize_name(name)
            scores = {}

            # Always try to get pLDDT from CIF
            cif_path = os.path.join(chai_dir, safe_name, "pred.model_idx_0.cif")
            plddt = extract_chai_confidence(cif_path)
            if plddt is not None:
                scores["plddt_mean"] = plddt

            # Try NPZ for pTM, iPTM, aggregate_score
            for npz_name in ["scores.model_idx_0.npz", "scores.model_idx_1.npz"]:
                npz_path = os.path.join(chai_dir, safe_name, npz_name)
                if os.path.exists(npz_path):
                    npz_scores = load_chai_scores(npz_path)
                    if npz_scores:
                        scores.update(npz_scores)
                        break

            if scores:
                chai_scores[safe_name] = scores

        mpnn_scores = {}
        for fname in sorted(os.listdir(mpnn_dir)):
            fpath = os.path.join(mpnn_dir, fname)
            if os.path.isdir(fpath):
                # Look for NPZ files first
                npz_files = [f for f in os.listdir(fpath) if f.endswith(".npz")]
                # Also check score_only subdirectory
                score_only_dir = os.path.join(fpath, "score_only")
                if os.path.isdir(score_only_dir):
                    npz_files.extend([os.path.join("score_only", f) for f in os.listdir(score_only_dir) if f.endswith(".npz")])

                if npz_files:
                    # Use first NPZ file
                    best = extract_mpnn_scores_from_npz(os.path.join(fpath, npz_files[0]))
                    if best is not None:
                        mpnn_scores[fname.replace("mpnn_", "")] = best
                else:
                    # Fallback to FASTA files
                    seqs_dir = os.path.join(fpath, "seqs")
                    if os.path.isdir(seqs_dir):
                        fas = [f for f in os.listdir(seqs_dir) if f.endswith(".fa")]
                        if fas:
                            s = extract_mpnn_scores(os.path.join(seqs_dir, fas[0]))
                            if s:
                                mpnn_scores[fname.replace("mpnn_", "")] = min(s.values())

        report_path = os.path.join(root, "chimeric_report.csv")
        ranked = write_full_report(report_path, chimeras, chimera_loops, query_seq,
                                    loops, chai_scores, mpnn_scores)
        print(f"  Done. Report: {report_path}")
        return

    # -----------------------------------------------------------------------
    # 1. Read inputs
    # -----------------------------------------------------------------------
    print("=" * 60)
    print("STEP 1: Reading input files")
    print("=" * 60)
    _, query_seq = read_fasta(args.fasta)
    print(f"  PH20M3 sequence length: {len(query_seq)}")

    msa_entries = parse_a3m(args.msa)
    print(f"  MSA entries: {len(msa_entries)}")

    # -----------------------------------------------------------------------
    # 2. Parse loops
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 2: Defining loops")
    print("=" * 60)
    loops = parse_loops(args.loops)
    for q_start, q_end in loops:
        print(f"  Loop {q_start}-{q_end}: {query_seq[q_start-1:q_end]}")

    # -----------------------------------------------------------------------
    # 3. Find human homologs
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 3: Finding Homo sapiens homologs (identity > {})".format(args.min_identity))
    print("=" * 60)
    # The first MSA entry is the query itself
    query_aln = msa_entries[0][1]
    homologs = find_human_homologs(msa_entries[1:], query_aln, args.min_identity)
    print(f"  Found {len(homologs)} human homolog sequences")
    for hdr, aln_seq, ident, label in homologs:
        print(f"    {label}: identity={ident:.3f}")

    if not homologs:
        print("  ERROR: No human homologs found. Cannot build chimeras.")
        sys.exit(1)

    # -----------------------------------------------------------------------
    # 4. Build chimeras (combinatorial)
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 4: Building chimeric sequences (combinatorial)")
    print("=" * 60)
    chimeras, chimera_loops = build_chimeras(query_seq, homologs, loops, args.max_chimeras)
    print(f"  Generated {len(chimeras)} unique chimeric sequences")

    if not chimeras:
        print("  ERROR: No chimeras generated. Check loop regions and MSA coverage.")
        sys.exit(1)

    # -----------------------------------------------------------------------
    # 5. Write FASTA files
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 5: Writing FASTA files for Chai-1")
    print("=" * 60)
    chai_fastas = []
    for name, seq in chimeras:
        safe_name = sanitize_name(name)
        fpath = os.path.join(fasta_dir, f"{safe_name}.fasta")
        write_fasta(seq, fpath, name)
        chai_fastas.append((name, safe_name, fpath))
        print(f"  {name}: {len(seq)} aa -> {fpath}")

    # Also write a combined FASTA with all chimeras
    combined_fasta = os.path.join(root, "all_chimeras.fasta")
    with open(combined_fasta, "w") as fh:
        for name, seq in chimeras:
            fh.write(f">protein|{name}\n{seq}\n")
    print(f"\n  Combined FASTA: {combined_fasta}")

    # Write a CSV index of chimeras
    index_csv = os.path.join(root, "chimeras_index.csv")
    loop_cols = []
    for q_start, q_end in loops:
        loop_cols.append(f"loop_{q_start}-{q_end}_donor")
        loop_cols.append(f"loop_{q_start}-{q_end}_seq")
    with open(index_csv, "w") as fh:
        fh.write("name,safe_name,sequence,length,fasta_path," + ",".join(loop_cols) + "\n")
        for name, seq in chimeras:
            safe_name = sanitize_name(name)
            loop_info = chimera_loops.get(name, [])
            loop_dict = {(s, e): (lbl, ds) for s, e, lbl, ds in loop_info}
            loop_cells = []
            for q_start, q_end in loops:
                if (q_start, q_end) in loop_dict:
                    lbl, ds = loop_dict[(q_start, q_end)]
                    loop_cells.append(lbl)
                    loop_cells.append(ds)
                else:
                    loop_cells.append("PH20M3_wt")
                    loop_cells.append(query_seq[q_start-1:q_end])
            fh.write(f"{name},{safe_name},{seq},{len(seq)},"
                     f"{os.path.join(fasta_dir, safe_name + '.fasta')},"
                     + ",".join(loop_cells) + "\n")
    print(f"  Chimera index: {index_csv}")

    # -----------------------------------------------------------------------
    # 6 & 7. Chai-1 prediction
    # -----------------------------------------------------------------------
    if not args.skip_chai1 and not args.score_only:
        print("\n" + "=" * 60)
        print("STEP 6: Running Chai-1 structure prediction")
        print("=" * 60)

        # Write a SLURM batch script for all jobs
        batch_script = os.path.join(root, "submit_all_chai1.sh")
        with open(batch_script, "w") as fh:
            fh.write("#!/bin/bash\n")
            fh.write(f"# Batch Chai-1 submission for {len(chai_fastas)} chimeras\n\n")
            for name, safe_name, fpath in chai_fastas:
                out_dir = os.path.join(chai_dir, safe_name)
                os.makedirs(out_dir, exist_ok=True)
                if args.submit_slurm:
                    chai_cmd = f"{CHAI1_RUN} {fpath} {out_dir}"
                    if args.msa_for_chai1:
                        chai_cmd += f" {os.path.abspath(args.msa_for_chai1)}"
                    fh.write(f"# {name}\n")
                    fh.write(f"{SLURM_SUBMIT} \"{chai_cmd}\" {args.ncpus} {args.slurm_partition}\n")
                else:
                    fh.write(f"# {name}\n")
                    cmd_str = f"{CHAI1_RUN} {fpath} {out_dir}"
                    if args.msa_for_chai1:
                        cmd_str += f" {os.path.abspath(args.msa_for_chai1)}"
                    fh.write(f"{cmd_str}\n")

        os.chmod(batch_script, 0o755)
        print(f"  Batch submission script: {batch_script}")

        if args.submit_slurm:
            # Submit each job individually (for better scheduling)
            for name, safe_name, fpath in chai_fastas:
                out_dir = os.path.join(chai_dir, safe_name)
                os.makedirs(out_dir, exist_ok=True)
                print(f"  Submitting Chai-1 for {name} ...")
                submit_chai1(fpath, out_dir,
                             msa_path=os.path.abspath(args.msa_for_chai1) if args.msa_for_chai1 else None,
                             submit_slurm=True,
                             slurm_partition=args.slurm_partition,
                             ncpus=args.ncpus)
        else:
            print("  To run Chai-1, execute the batch script or use --submit-slurm")
            print(f"    bash {batch_script}")

    # -----------------------------------------------------------------------
    # 7. Convert CIF → PDB
    # -----------------------------------------------------------------------
    cif_available = []
    if not args.score_only and not args.skip_chai1 and not args.report_only:
        print("\n" + "=" * 60)
        print("STEP 7: Converting CIF → PDB (model_idx_0)")
        print("=" * 60)
        make_cif2pdb_script()

        for name, _ in chimeras:
            safe_name = sanitize_name(name)
            cif_path = os.path.join(chai_dir, safe_name, "pred.model_idx_0.cif")
            pdb_path = os.path.join(pdb_dir, f"{safe_name}.pdb")
            if os.path.exists(cif_path):
                ok = convert_cif_to_pdb(cif_path, pdb_path)
                if ok:
                    cif_available.append((safe_name, cif_path))
                    print(f"  {name}: CIF -> PDB ok")
            else:
                print(f"  {name}: CIF not found (prediction may not be complete yet)")
                print(f"    Expected: {cif_path}")

    # -----------------------------------------------------------------------
    # 8. Extract Chai-1 scores from NPZ / CIF
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 8: Extracting Chai-1 confidence scores")
    print("=" * 60)

    chai_scores = {}  # safe_name -> {plddt_mean, ptm, iptm, ...}
    for name, _ in chimeras:
        safe_name = sanitize_name(name)
        scores = {}

        # Always try to get pLDDT from CIF
        cif_path = os.path.join(chai_dir, safe_name, "pred.model_idx_0.cif")
        plddt = extract_chai_confidence(cif_path)
        if plddt is not None:
            scores["plddt_mean"] = plddt

        # Try NPZ for pTM, iPTM, aggregate_score
        for npz_name in ["scores.model_idx_0.npz", "scores.model_idx_1.npz"]:
            npz_path = os.path.join(chai_dir, safe_name, npz_name)
            if os.path.exists(npz_path):
                npz_scores = load_chai_scores(npz_path)
                if npz_scores:
                    scores.update(npz_scores)
                    break

        if scores:
            chai_scores[safe_name] = scores

    n_scored = len(chai_scores)
    n_plddt = sum(1 for s in chai_scores if "plddt_mean" in chai_scores[s])
    n_ptm = sum(1 for s in chai_scores if "ptm" in chai_scores[s])
    print(f"  Chai-1 scores extracted: {n_scored}/{len(chimeras)} (pLDDT: {n_plddt}, pTM: {n_ptm})")

    # -----------------------------------------------------------------------
    # 9. Score with ProteinMPNN
    # -----------------------------------------------------------------------
    mpnn_scores = {}  # safe_name -> best MPNN score
    if not args.skip_mpnn:
        print("\n" + "=" * 60)
        print("STEP 9: Scoring PDBs with ProteinMPNN")
        print("=" * 60)

        # Find PDBs: from pdb_dir, or convert CIFs on the fly
        pdb_files = sorted(os.listdir(pdb_dir)) if os.path.isdir(pdb_dir) else []
        if not pdb_files:
            for name, _ in chimeras:
                safe_name = sanitize_name(name)
                cif_path = os.path.join(chai_dir, safe_name, "pred.model_idx_0.cif")
                pdb_path = os.path.join(pdb_dir, f"{safe_name}.pdb")
                if os.path.exists(cif_path) and not os.path.exists(pdb_path):
                    convert_cif_to_pdb(cif_path, pdb_path)
            pdb_files = sorted(os.listdir(pdb_dir)) if os.path.isdir(pdb_dir) else []

        if not pdb_files:
            print("  No PDB files found for MPNN scoring.")
        else:
            # Create mapping from safe_name to sequence
            seq_map = {sanitize_name(name): seq for name, seq in chimeras}

            for pdb_file in pdb_files:
                if not pdb_file.endswith(".pdb"):
                    continue
                pdb_path = os.path.join(pdb_dir, pdb_file)
                seq_name = pdb_file.replace(".pdb", "")
                print(f"  Scoring {seq_name} ...")
                chimera_seq = seq_map.get(seq_name)
                output_path = score_pdb_with_proteinmpnn(pdb_path, mpnn_dir, seq_name, chimera_seq)
                if output_path:
                    if output_path.endswith('.npz'):
                        # NPZ file format
                        best = extract_mpnn_scores_from_npz(output_path)
                        if best is not None:
                            mpnn_scores[seq_name] = best
                            print(f"    Best MPNN score: {best:.4f}")
                    else:
                        # FASTA file format (legacy)
                        scores = extract_mpnn_scores(output_path)
                        if scores:
                            best = min(scores.values())
                            mpnn_scores[seq_name] = best
                            print(f"    Best MPNN score: {best:.4f}")

    # -----------------------------------------------------------------------
    # 10. Write comprehensive CSV report
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 10: Writing comprehensive CSV report")
    print("=" * 60)
    report_path = os.path.join(root, "chimeric_report.csv")
    ranked = write_full_report(report_path, chimeras, chimera_loops, query_seq,
                                loops, chai_scores, mpnn_scores)

    # Print top-10 summary
    print("\n  Top 10 chimeras by combined score:")
    print(f"  {'Rank':<6} {'Name':<50} {'pLDDT':<8} {'pTM':<8} {'iPTM':<8} {'MPNN':<8}")
    print("  " + "-" * 88)
    for row in ranked[:10]:
        print(f"  {row['rank_combined']:<6} {row['name'][:48]:<50} "
              f"{row['chai1_plddt']:<8} {row['chai1_ptm']:<8} "
              f"{row['chai1_iptm']:<8} {row['mpnn_best_score']:<8}")

    # Also write the old-style score summary if MPNN scores exist
    if mpnn_scores:
        summary_path = os.path.join(root, "mpnn_scores_summary.csv")
        write_score_summary(mpnn_scores, summary_path)

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)
    print(f"  Output directory: {root}")
    print(f"  Chimeric sequences: {len(chimeras)}")
    print(f"  Chai-1 scored:     {n_scored}/{len(chimeras)}")
    print(f"  MPNN scored:       {len(mpnn_scores)}/{len(chimeras)}")
    print(f"  FASTA files:       {fasta_dir}")
    print(f"  Chai-1 results:    {chai_dir}")
    print(f"  PDB files:         {pdb_dir}")
    print(f"  MPNN scores:       {mpnn_dir}")
    print(f"  Full report:       {report_path}")
    print(f"  Combined FASTA:    {combined_fasta}")
    print(f"  Chimera index:     {index_csv}")
    print()


if __name__ == "__main__":
    main()
