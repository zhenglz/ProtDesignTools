#!/usr/bin/env python3
"""
RFDiffusion3 Design Tool
========================

A tool for batch RFDiffusion3 protein design. Accepts PDB files and design region
specifications, generates JSON configurations, and submits design jobs.

All outputs (designs, scores, top sequences) are placed directly in the output
directory for easy inspection.

Design region syntax:
  <chain><start>-<end>:<design_len_min>-<design_len_max>[;<next_region>...]

  Replace existing residues — keep the surrounding context, design a new segment
  of the specified length in place of the original residues.

  Examples:
    --design-regions "A5-20:10-25"
        Replace chain A residues 5-20 (1-based) with a new 10-25 aa segment

    --design-regions "A90-95:5-5;B167-180:10-15"
        Replace A90-95 with 5 aa, replace B167-180 with 10-15 aa

    --design-regions "A0-0:5-10"
        Add 5-10 aa N-terminal extension on chain A

Usage:
  python tools/rfdiffusion3_tool.py --pdb ./input.pdb --output ./designs \\
    --design-regions "A5-20:10-25" --num-designs 10 --top-n 3

  python tools/rfdiffusion3_tool.py --pdb-dir ./input_pdbs --output ./designs \\
    --design-regions "A90-95:5-5;B167-180:10-15" --num-designs 5

  python tools/rfdiffusion3_tool.py --pdb ./input.pdb --output ./designs \\
    --design-regions "A0-0:5-10" --local
"""

import argparse
import csv
import glob
import gzip
import json
import os
import re
import shutil
import subprocess as sp
import sys
import time
from collections import OrderedDict


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
RFD3_BIN = "/data_test/home/guoliangzhu/miniconda3/envs/rfd3/bin/rfd3"
RFD3_EXAMPLE_DIR = "/data_test/home/guoliangzhu/bioapp/rfdiffusion3/example"
SLURM_SUBMIT = "/data_test/home/lzzheng/bin/submit_slurm_gpu.sh"


# ---------------------------------------------------------------------------
# 3-letter to 1-letter amino acid code
# ---------------------------------------------------------------------------
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

# Standard DNA residues (BioPython convention)
DNA_RESIDUES = {'DA', 'DT', 'DC', 'DG', 'A', 'T', 'C', 'G'}
# Standard RNA residues
RNA_RESIDUES = {'A', 'U', 'C', 'G', 'RA', 'RU', 'RC', 'RG'}
# Combined set of everything that is NOT protein
NON_PROTEIN = DNA_RESIDUES | RNA_RESIDUES | {'ADE', 'THY', 'CYT', 'GUA', 'URA'}


def parse_args():
    parser = argparse.ArgumentParser(
        description="RFDiffusion3 Design Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Input sources (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--pdb-dir", help="Directory containing input PDB files")
    input_group.add_argument("--pdb", help="Single input PDB file")
    input_group.add_argument("--json", help="JSON config file for RFDiffusion3")

    # Output
    parser.add_argument("--output", "-o", default="./rfd3_designs",
                        help="Output directory (default: ./rfd3_designs)")

    # Design regions
    parser.add_argument("--design-regions", default=None,
                        help="Design region specification. "
                             "Format: <chain><start>-<end>:<len_min>-<len_max>"
                             "[;<chain><start>-<end>:<len_min>-<len_max>...]")

    # Explicit contig (overrides --design-regions)
    parser.add_argument("--contig", default=None,
                        help="Explicit contig string (overrides --design-regions). "
                             "E.g. 'A1-4/10-25,A21-150'")
    parser.add_argument("--fixed-atoms", default=None,
                        help="Fixed atoms JSON string or 'none'. Overrides auto-fixing.")
    parser.add_argument("--num-designs", type=int, default=1,
                        help="Number of designs per input (default: 1).")

    # Top-N selection
    parser.add_argument("--top-n", type=int, default=0,
                        help="Number of top designs to include in the summary FASTA "
                             "(0 = all). Ranked by a combined score.")

    # Execution mode
    parser.add_argument("--local", action="store_true",
                        help="Run locally instead of SLURM")
    parser.add_argument("--slurm-partition", default="4090",
                        help="SLURM partition (default: 4090)")
    parser.add_argument("--ncpus", type=int, default=4,
                        help="CPUs per SLURM job (default: 4)")

    # Concurrency
    parser.add_argument("--max-jobs", type=int, default=4,
                        help="Maximum concurrent SLURM jobs (default: 4)")

    # Processing
    parser.add_argument("--no-wait", action="store_true",
                        help="Submit jobs and exit without waiting")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip PDBs/JSONs with existing outputs")
    parser.add_argument("--extract-only", action="store_true",
                        help="Re-process existing outputs (score extraction, FASTA)")

    return parser.parse_args()


# ---------------------------------------------------------------------------
# PDB helpers
# ---------------------------------------------------------------------------

def find_pdb_files(pdb_dir):
    """Find all .pdb files in a directory recursively."""
    return glob.glob(os.path.join(pdb_dir, "**", "*.pdb"), recursive=True)


def guess_chain_length(pdb_file, chain_id):
    """Return the last residue number (1-based) for a chain in the PDB."""
    max_res = 0
    try:
        with open(pdb_file) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    cid = line[21].strip()
                    if cid == chain_id:
                        try:
                            resi = int(line[22:26].strip())
                            if resi > max_res:
                                max_res = resi
                        except ValueError:
                            pass
        return max_res
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Design regions → contig
# ---------------------------------------------------------------------------

def parse_design_regions(spec):
    """Parse 'A5-20:10-25;B167-180:10-15' into structured specs."""
    if not spec:
        return None
    regions = []
    pattern = re.compile(r'^([A-Za-z])(\d+)-(\d+):(\d+)-(\d+)$')
    for part in spec.split(';'):
        part = part.strip()
        m = pattern.match(part)
        if not m:
            raise ValueError(
                f"Cannot parse design region: '{part}'. "
                f"Expected: <chain><start>-<end>:<len_min>-<len_max>"
            )
        regions.append({
            'chain': m.group(1),
            'start': int(m.group(2)),
            'end': int(m.group(3)),
            'len_min': int(m.group(4)),
            'len_max': int(m.group(5)),
        })
    return regions


def build_contig(regions, pdb_file=None):
    """Build an RFDiffusion3 contig string from parsed design region specs.

    For each region (start > 0): keep context-before, design /Lmin-Lmax.
    For N-terminal (start == 0): design /Lmin-Lmax, keep entire chain.
    Trailing context added after all regions on each chain.
    """
    by_chain = OrderedDict()
    for r in regions:
        by_chain.setdefault(r['chain'], []).append(r)
    for rlist in by_chain.values():
        rlist.sort(key=lambda x: x['start'])

    parts = []
    for chain, rlist in by_chain.items():
        chain_len = guess_chain_length(pdb_file, chain) if pdb_file else 999

        prev_end = 0  # last residue *kept* before a design region

        for r in rlist:
            s, e = r['start'], r['end']

            if s == 0:
                # N-terminal addition: design then the whole chain
                parts.append(f"/{r['len_min']}-{r['len_max']}")
                parts.append(f"{chain}1-{chain_len}")
                prev_end = chain_len
            else:
                # Context before (if there is a gap since prev_end)
                if s > prev_end + 1:
                    parts.append(f"{chain}{prev_end + 1}-{s - 1}")
                # Design segment replacing s..e
                parts.append(f"/{r['len_min']}-{r['len_max']}")
                prev_end = e

        # Trailing context after last region on this chain
        if prev_end < chain_len:
            parts.append(f"{chain}{prev_end + 1}-{chain_len}")

    contig = ",".join(parts).replace(",/", "/")
    contig = re.sub(r'[A-Za-z]1-0/', '', contig)
    return contig


# ---------------------------------------------------------------------------
# JSON preparation
# ---------------------------------------------------------------------------

def prepare_json(pdb_file, output_prefix, design_regions=None, contig_str=None,
                 fixed_atoms=None, fixed_all_context=True):
    """Create an RFDiffusion3 JSON config and return its path."""
    out_dir = os.path.dirname(output_prefix)
    os.makedirs(out_dir, exist_ok=True)

    # Build contig from design regions
    if design_regions and not contig_str:
        contig_str = build_contig(design_regions, pdb_file)

    # Fixed atoms — by default fix everything except the designed regions
    fixed_atoms_dict = {}
    if fixed_atoms:
        if isinstance(fixed_atoms, str):
            try:
                fixed_atoms_dict = json.loads(fixed_atoms)
            except json.JSONDecodeError:
                fixed_atoms_dict = fixed_atoms if fixed_atoms.lower() != "none" else {}
        else:
            fixed_atoms_dict = fixed_atoms
    elif fixed_all_context and design_regions:
        for r in design_regions:
            ch, s, e = r['chain'], r['start'], r['end']
            chain_len = guess_chain_length(pdb_file, ch) or 9999
            if s > 0:
                if s > 1:
                    fixed_atoms_dict[f"{ch}1-{s - 1}"] = "ALL"
                if e < chain_len:
                    fixed_atoms_dict[f"{ch}{e + 1}-{chain_len}"] = "ALL"
            else:
                fixed_atoms_dict[f"{ch}1-{chain_len}"] = "ALL"

    json_data = {
        "protein_design": {
            "input": os.path.abspath(pdb_file),
            "contig": contig_str or "",
            "length": "",
            "select_fixed_atoms": fixed_atoms_dict,
            "is_non_loopy": True,
            "dialect": 2,
        }
    }

    json_path = f"{output_prefix}.json"
    with open(json_path, 'w') as f:
        json.dump(json_data, f, indent=2)
    return json_path


# ---------------------------------------------------------------------------
# Job submission
# ---------------------------------------------------------------------------

def submit_slurm(cmd, ncpus, partition):
    """Submit a command via SLURM GPU submission script."""
    if not os.path.exists(SLURM_SUBMIT):
        print(f"Error: SLURM submit script not found: {SLURM_SUBMIT}")
        return None
    try:
        result = sp.run(
            [SLURM_SUBMIT, cmd, str(ncpus), partition],
            capture_output=True, text=True, check=True,
        )
        out = result.stdout.strip()
        for line in out.split('\n'):
            if "Submitted batch job" in line:
                return line.split()[-1]
        if out.isdigit():
            return out
        print(f"Warning: Could not parse job ID from: {out}")
        return None
    except sp.CalledProcessError as e:
        print(f"SLURM submission failed: {e.stderr}")
        return None


def run_local(cmd, job_name, log_dir):
    """Run a command locally in background, return PID."""
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"{job_name}.log")
    try:
        f = open(log_file, 'w')
        p = sp.Popen(cmd, shell=True, stdout=f, stderr=sp.STDOUT)
        f.close()
        return str(p.pid)
    except Exception as e:
        print(f"Local run failed for {job_name}: {e}")
        return None


def check_job_status_slurm(job_id):
    """Check SLURM job status via squeue."""
    if not job_id:
        return "COMPLETED"
    try:
        r = sp.run(['squeue', '-j', str(job_id), '-h', '-o', '%T'],
                   capture_output=True, text=True)
        if r.returncode == 0:
            state = r.stdout.strip()
            return state if state else "COMPLETED"
        return "COMPLETED"
    except Exception:
        return "UNKNOWN"


def check_job_status_local(pid):
    """Check local process status."""
    if not pid:
        return "COMPLETED"
    try:
        r = sp.run(['ps', '-p', str(pid), '-o', 'state', '--noheaders'],
                   capture_output=True, text=True)
        return "RUNNING" if r.returncode == 0 else "COMPLETED"
    except Exception:
        return "UNKNOWN"


# ---------------------------------------------------------------------------
# Output processing
# ---------------------------------------------------------------------------

def guess_chain_type(residues):
    """Classify a chain as 'protein', 'dna', 'rna', or 'ligand'.

    `residues` is a list of BioPython residue objects.
    """
    resnames = set(r.resname for r in residues if r.id[0] == ' ')
    if not resnames:
        return 'ligand'

    # Count how many residues match each type
    n_protein = sum(1 for rn in resnames if rn in THREE_TO_ONE)
    n_dna = sum(1 for rn in resnames if rn in DNA_RESIDUES)
    n_ligand = sum(1 for rn in resnames
                   if rn not in THREE_TO_ONE and rn not in NON_PROTEIN)

    if n_dna > 0 and n_dna >= n_protein and n_dna >= n_ligand:
        return 'dna'
    if n_protein >= n_dna and n_protein >= n_ligand:
        return 'protein'
    if n_dna > 0:
        return 'dna'
    return 'ligand'


def one_letter_seq(residues):
    """Convert a list of BioPython residues to a 1-letter sequence.

    Protein residues → 1-letter amino acid.
    DNA residues → standard letters (A, T, C, G).
    RNA residues → standard letters (A, U, C, G).
    Others → 'X'.
    """
    seq = []
    for r in residues:
        if r.id[0] != ' ':
            continue
        rn = r.resname
        if rn in THREE_TO_ONE:
            seq.append(THREE_TO_ONE[rn])
        elif rn in DNA_RESIDUES:
            if rn.startswith('D'):
                seq.append(rn[1])  # DA → A, DT → T
            else:
                seq.append(rn[0])  # A → A
        elif rn in RNA_RESIDUES:
            seq.append(rn[0])
        else:
            seq.append('X')
    return ''.join(seq)


def cif_sequence(cif_file):
    """Extract sequences from a CIF file, classifying each chain.

    Returns a list of (chain_id, chain_type, sequence_string) tuples.
    chain_type is one of 'protein', 'dna', 'rna', 'ligand'.
    """
    try:
        from Bio.PDB import MMCIFParser
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("model", cif_file)
        chains = []
        for model in structure:
            for chain in model:
                residues = [r for r in chain if r.id[0] == ' ']
                if not residues:
                    continue
                chain_type = guess_chain_type(residues)
                seq = one_letter_seq(residues)
                if seq:
                    chains.append((chain.id, chain_type, seq))
        return chains
    except ImportError:
        return None
    except Exception:
        return None


def load_design_scores(json_path):
    """Load metrics from an RFDiffusion3 output JSON.

    Returns a dict with keys: max_ca_deviation, n_chainbreaks, n_clashing,
    loop_fraction, helix_fraction, sheet_fraction, radius_of_gyration,
    num_residues, alanine_content, glycine_content.
    """
    try:
        with open(json_path) as f:
            data = json.load(f)
        return data.get("metrics", data)
    except (json.JSONDecodeError, IOError, KeyError):
        return {}


def process_outputs(task_dir, output_dir):
    """Process a single task's outputs.

    Steps:
      1. Decompress .cif.gz → .cif
      2. Convert .cif → .pdb (via BioPython)
      3. Extract sequences from CIF
      4. Load design scores from companion JSON

    All files remain in their original task directory. Returns a list of
    dicts with design metadata.
    """
    gz_files = sorted(glob.glob(os.path.join(task_dir, "*.cif.gz")))
    designs = []

    for gz_path in gz_files:
        base = os.path.splitext(os.path.basename(gz_path))[0]  # e.g. demo_0_model_0.cif
        if base.endswith('.cif'):
            base = base[:-4]

        cif_path = os.path.join(task_dir, f"{base}.cif")
        pdb_path = os.path.join(task_dir, f"{base}.pdb")
        json_path = os.path.join(task_dir, f"{base}.json")

        # Decompress CIF
        if not os.path.exists(cif_path) and os.path.exists(gz_path):
            try:
                with gzip.open(gz_path, 'rb') as f_in:
                    with open(cif_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
            except Exception as e:
                print(f"    Error decompressing {gz_path}: {e}")

        # Convert CIF → PDB
        if os.path.exists(cif_path) and not os.path.exists(pdb_path):
            try:
                from Bio.PDB import MMCIFParser, PDBIO
                parser = MMCIFParser(QUIET=True)
                structure = parser.get_structure("model", cif_path)
                io = PDBIO()
                io.set_structure(structure)
                io.save(pdb_path)
            except ImportError:
                print("    BioPython not available, skipping CIF→PDB conversion")
            except Exception as e:
                print(f"    Error converting {cif_path}: {e}")

        # Extract sequences
        chains = None
        if os.path.exists(cif_path):
            chains = cif_sequence(cif_path)

        # Load scores
        scores = {}
        if os.path.exists(json_path):
            scores = load_design_scores(json_path)

        entry = {
            "design_name": base,
            "task_dir": task_dir,
            "cif_file": cif_path if os.path.exists(cif_path) else None,
            "pdb_file": pdb_path if os.path.exists(pdb_path) else None,
            "json_file": json_path if os.path.exists(json_path) else None,
            "chains": chains,
            "scores": scores,
            # Composite "badness" score — lower is better
            "max_ca_deviation": scores.get("max_ca_deviation", 0),
            "n_chainbreaks": scores.get("n_chainbreaks", 0),
            "n_clashing": scores.get("n_clashing.interresidue_clashes_w_sidechain", 0)
                          or scores.get("n_clashing", 0),
            "loop_fraction": scores.get("loop_fraction", 0),
            "helix_fraction": scores.get("helix_fraction", 0),
            "sheet_fraction": scores.get("sheet_fraction", 0),
            "num_residues": scores.get("num_residues", 0),
            "radius_of_gyration": scores.get("radius_of_gyration", 0),
            "alanine_content": scores.get("alanine_content", 0),
            "glycine_content": scores.get("glycine_content", 0),
        }
        designs.append(entry)

    return designs


def composite_rank(design):
    """Compute a ranking score: lower = better design.

    Penalises: chainbreaks, clashes, extreme RG, high alanine/glycine.
    """
    score = (
        design["n_chainbreaks"] * 100 +
        design["n_clashing"] * 10 +
        abs(design["radius_of_gyration"] - 12) * 0.1 +
        design["alanine_content"] * 10 +
        abs(design["glycine_content"] - 0.05) * 10 +
        design["max_ca_deviation"]
    )
    return score


def write_fasta(designs, output_dir, top_n=0):
    """Write per-design multi-chain FASTA files and a combined FASTA.

    Each top design gets its own FASTA file (all chains together) in a
    'top_designs_fastas/' subdirectory. A combined FASTA with all selected
    designs is also written at the output root.

    Headers use the Chai-1 convention:
      >protein|design_01_A      for protein chains
      >DNA|design_01_B          for DNA chains
      >ligand|design_01_C       for ligand / other chains

    If top_n > 0, only the top N designs (by composite score) are included.
    """
    if not designs:
        return

    ranked = sorted(designs, key=composite_rank)
    if top_n > 0:
        ranked = ranked[:top_n]

    fasta_dir = os.path.join(output_dir, "top_designs_fastas")
    os.makedirs(fasta_dir, exist_ok=True)

    for rank, d in enumerate(ranked, 1):
        if not d["chains"]:
            continue
        design_tag = f"design_{rank:02d}"
        fasta_path = os.path.join(fasta_dir, f"{d['design_name']}.fasta")
        with open(fasta_path, 'w') as f:
            for cid, chain_type, seq in d["chains"]:
                header = (
                    f"{chain_type_prefix(chain_type)}|{design_tag}_{cid} "
                    f"ca_dev={d['max_ca_deviation']:.2f} "
                    f"breaks={d['n_chainbreaks']} "
                    f"clashes={d['n_clashing']} "
                    f"rg={d['radius_of_gyration']:.1f}"
                )
                f.write(f"{header}\n{seq}\n")

    # Combined FASTA at the root
    combined_fasta = os.path.join(output_dir, "top_designs.fasta")
    with open(combined_fasta, 'w') as f:
        for rank, d in enumerate(ranked, 1):
            if not d["chains"]:
                continue
            design_tag = f"design_{rank:02d}"
            for cid, chain_type, seq in d["chains"]:
                header = (
                    f"{chain_type_prefix(chain_type)}|{design_tag}_{cid} "
                    f"ca_dev={d['max_ca_deviation']:.2f} "
                    f"breaks={d['n_chainbreaks']} "
                    f"clashes={d['n_clashing']} "
                    f"rg={d['radius_of_gyration']:.1f}"
                )
                f.write(f"{header}\n{seq}\n")

    n = len(ranked)
    print(f"  Wrote {n} {'top ' if top_n else ''}design{'s' if n != 1 else ''} to {fasta_dir}/")
    print(f"  Combined FASTA: {combined_fasta}")


def chain_type_prefix(chain_type):
    """Map chain type to FASTA header prefix (Chai-1 convention)."""
    return {
        'protein': 'protein',
        'dna': 'DNA',
        'rna': 'RNA',
    }.get(chain_type, 'ligand')


def write_summary_csv(all_designs, csv_path):
    """Write a CSV summary of all designs with scores."""
    fields = [
        "design_name", "task_dir",
        "max_ca_deviation", "n_chainbreaks", "n_clashing",
        "loop_fraction", "helix_fraction", "sheet_fraction",
        "num_residues", "radius_of_gyration",
        "alanine_content", "glycine_content",
        "cif_file", "pdb_file",
    ]
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        for d in all_designs:
            row = {k: d.get(k, d["scores"].get(k, "")) for k in fields}
            row["design_name"] = d["design_name"]
            row["task_dir"] = d["task_dir"]
            row["cif_file"] = d["cif_file"]
            row["pdb_file"] = d["pdb_file"]
            w.writerow(row)
    print(f"  Wrote summary CSV: {csv_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # ---- Phase 1: Gather inputs & prepare JSONs ----
    json_files = []

    if args.json:
        json_files = [os.path.abspath(args.json)]
        print(f"Using JSON config: {json_files[0]}")
    else:
        if args.pdb_dir:
            pdbs = find_pdb_files(args.pdb_dir)
            if not pdbs:
                pdbs = glob.glob(os.path.join(args.pdb_dir, "*.pdb"))
            if not pdbs:
                print(f"No PDB files found in {args.pdb_dir}")
                sys.exit(1)
        elif args.pdb:
            pdbs = [os.path.abspath(args.pdb)]
        else:
            print("No input specified")
            sys.exit(1)

        print(f"Found {len(pdbs)} PDB file(s).")

        design_regions = None
        if args.design_regions:
            design_regions = parse_design_regions(args.design_regions)
            region_str = "; ".join(
                f"{r['chain']}{r['start']}-{r['end']}:{r['len_min']}-{r['len_max']}"
                for r in design_regions
            )
            print(f"Design regions: {region_str}")

        # Parse fixed atoms
        fixed_atoms = None
        if args.fixed_atoms and args.fixed_atoms.lower() != "none":
            try:
                fixed_atoms = json.loads(args.fixed_atoms)
            except json.JSONDecodeError:
                print("  Warning: Could not parse --fixed-atoms JSON, ignoring")

        print("\nPreparing RFDiffusion3 JSON configurations...")
        for pdb in pdbs:
            pdb_name = os.path.splitext(os.path.basename(pdb))[0]
            # All outputs for this PDB go directly in the output directory
            # with the PDB name as prefix
            prefix = os.path.join(args.output, pdb_name)
            os.makedirs(args.output, exist_ok=True)

            json_path = prepare_json(
                pdb, prefix,
                design_regions=design_regions,
                contig_str=args.contig,
                fixed_atoms=fixed_atoms,
                fixed_all_context=True,
            )
            json_files.append(json_path)
            print(f"  {pdb_name}: {json_path}")

    if not json_files:
        print("No JSON configurations to process.")
        sys.exit(1)

    print(f"Total JSON configurations: {len(json_files)}")

    # ---- Phase 2: Filter for skip-existing ----
    if args.skip_existing and not args.extract_only:
        filtered = []
        for jf in json_files:
            base = os.path.splitext(os.path.basename(jf))[0]
            task_dir = args.output
            cif_gz = glob.glob(os.path.join(task_dir, f"{base}_model_*.cif.gz"))
            if len(cif_gz) >= args.num_designs:
                print(f"  Skipping {base}: outputs already exist")
                continue
            filtered.append(jf)
        json_files_submit = filtered
        json_files_process = json_files
    elif args.extract_only:
        json_files_submit = []
        json_files_process = json_files
    else:
        json_files_submit = json_files
        json_files_process = json_files

    # ---- Phase 3: Submit jobs ----
    log_dir = os.path.join(args.output, "logs")
    jobs = []  # [(job_id, name, task_dir, json_file)]

    if json_files_submit:
        print(f"\nSubmitting {len(json_files_submit)} RFDiffusion3 job(s)...")
        for jf in json_files_submit:
            base_name = os.path.splitext(os.path.basename(jf))[0]
            abs_json = os.path.abspath(jf)
            abs_out = os.path.abspath(args.output)

            cmd = f"cd {RFD3_EXAMPLE_DIR} && {RFD3_BIN} design out_dir={abs_out} inputs={abs_json}"

            if args.local:
                job_id = run_local(cmd, base_name, log_dir)
            else:
                job_id = submit_slurm(cmd, args.ncpus, args.slurm_partition)

            if job_id:
                print(f"  {base_name}: submitted (ID: {job_id})")
            else:
                print(f"  {base_name}: submission FAILED")

            jobs.append((job_id, base_name, args.output, jf))
    else:
        for jf in json_files_process:
            base_name = os.path.splitext(os.path.basename(jf))[0]
            jobs.append((None, base_name, args.output, jf))

    if not jobs:
        print("Nothing to process.")
        return

    # ---- Phase 4: Wait for completion ----
    if args.no_wait:
        print("\nJobs submitted. Re-run without --no-wait (or with --extract-only) "
              "to collect results.")
        return

    submitted = [(jid, name, out_dir, jf) for jid, name, out_dir, jf in jobs if jid]
    if submitted:
        print(f"\nWaiting for {len(submitted)} job(s) to complete...")
        pending = {jid: (name, out_dir, jf) for jid, name, out_dir, jf in submitted}
        total = len(pending)
        completed_names = set()

        while pending:
            done = []
            for jid, (name, out_dir, jf) in pending.items():
                status = check_job_status_local(jid) if args.local else check_job_status_slurm(jid)
                if status in ("COMPLETED", "CD", "COMPLETING", "CG"):
                    done.append(jid)
                    completed_names.add(name)
                elif status in ("FAILED", "F", "CANCELLED", "CA", "TIMEOUT", "TO",
                                "NODE_FAIL", "NF", "BOOT_FAIL", "BF", "OUT_OF_MEMORY", "OOM"):
                    done.append(jid)
                    print(f"  {name}: FAILED ({status})")

            for jid in done:
                del pending[jid]

            if pending:
                ids_shown = list(pending.keys())[:3]
                done_count = total - len(pending)
                print(f"[{time.strftime('%H:%M:%S')}] {done_count}/{total} done, "
                      f"{len(pending)} pending (IDs: {', '.join(ids_shown)}{'...' if len(pending) > 3 else ''})")
                time.sleep(30)

        print(f"All {total} job(s) completed.")

    # ---- Phase 5: Process outputs ----
    print("\nProcessing outputs and collecting scores...")
    all_designs = []
    for jid, base_name, out_dir, jf in jobs:
        print(f"  {base_name}:")
        designs = process_outputs(out_dir, args.output)
        if designs:
            for d in designs:
                rank = composite_rank(d)
                d["_rank_score"] = rank
                print(f"    {d['design_name']}: "
                      f"CA_dev={d['max_ca_deviation']:.2f} "
                      f"breaks={d['n_chainbreaks']} "
                      f"clashes={d['n_clashing']} "
                      f"RG={d['radius_of_gyration']:.1f} "
                      f"n_res={d['num_residues']} "
                      f"score={rank:.2f}")
        else:
            print(f"    No outputs found")
        all_designs.extend(designs)

    # ---- Phase 6: Write summary CSV ----
    csv_path = os.path.join(args.output, "rfd3_scores.csv")
    write_summary_csv(all_designs, csv_path)

    # ---- Phase 7: Write per-design multi-chain FASTA files ----
    if all_designs:
        write_fasta(all_designs, args.output, top_n=args.top_n)
    else:
        print("No designs collected; skipping FASTA generation.")

    print("\nRFDiffusion3 design tool completed.")


if __name__ == "__main__":
    main()
