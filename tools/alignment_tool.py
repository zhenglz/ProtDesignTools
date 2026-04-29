#!/usr/bin/env python3
"""
Protein Structure Alignment Tool
=================================

Align two protein structures using TM-align.  Accepts PDB or mmCIF input
for both query and reference — CIF files are automatically converted to
PDB before alignment.

Outputs:
  1. TM-score (normalized by reference length), RMSD, and sequence identity
  2. Per-residue alignment mapping (1-indexed residue pairs)
  3. Aligned query structure (superimposed onto reference) as a PDB file
  4. All scores written to a JSON report

Tools and paths are defined in data/config.json under "alignment".

Usage:
  python tools/alignment_tool.py --query query.pdb --reference ref.pdb \\
    --output ./align_results

  python tools/alignment_tool.py --query query.cif --reference ref.pdb \\
    --output ./align_results

  python tools/alignment_tool.py --query query.pdb --reference ref.cif \\
    --output ./align_results \\
    --skip-existing
"""

import argparse
import json
import os
import re
import shutil
import subprocess as sp
import sys
import uuid

from _config import load_config


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run_command(cmd, desc=None):
    """Run a shell command, print status, and exit on failure."""
    tag = desc or cmd
    print(f"[RUN] {tag}")
    ret = sp.run(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                 universal_newlines=True)
    if ret.returncode != 0:
        print(f"[ERROR] Command failed (exit {ret.returncode}): {cmd}",
              file=sys.stderr)
        if ret.stderr:
            print(ret.stderr, file=sys.stderr)
        sys.exit(1)
    return ret


def get_file_ext(fpath):
    """Return lowercase extension without dot."""
    _, ext = os.path.splitext(fpath)
    return ext.lstrip(".").lower()


def is_cif(fpath):
    """Return True if the file is an mmCIF format."""
    return get_file_ext(fpath) in ("cif", "mmcif")


# ---------------------------------------------------------------------------
# CIF → PDB conversion wrapper
# ---------------------------------------------------------------------------

def cif_to_pdb(cif_path, pdb_path, config):
    """Convert a CIF file to PDB using cif2pdb.py."""
    python_exe = config["alignment"]["python_exe"]
    script = config["alignment"]["cif2pdb_script"]
    cmd = f"{python_exe} {script} {cif_path} {pdb_path}"
    run_command(cmd, "cif_to_pdb")
    if not os.path.exists(pdb_path):
        print(f"[ERROR] CIF→PDB conversion failed for {cif_path}",
              file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# TMalign output parser
# ---------------------------------------------------------------------------

def parse_tmalign_stdout(text):
    """Parse TM-align stdout and return a dict of scores.

    Extracts:
        TM-score (normalized by chain_1 / chain_2)
        RMSD
        aligned_length
        seq_id (n_identical / n_aligned)
        residue_pairs  (list of (res_index_in_query, res_index_in_ref, aa_q, aa_r, distance_flag))

    The residue-level alignment lines look like::

        AAGCTT
        ::::..
        AAGCTT

    where the middle line encodes spatial proximity (``:`` ≤ 5.0 Å, ``.`` > 5.0 Å).
    """
    result = {}

    # --- scores ---
    m = re.search(r"Aligned length=\s*(\d+),\s*RMSD=\s*([\d.]+),"
                  r"\s*Seq_ID=n_identical/n_aligned=\s*([\d.]+)",
                  text)
    if m:
        result["aligned_length"] = int(m.group(1))
        result["rmsd"] = float(m.group(2))
        result["seq_id"] = float(m.group(3))

    m = re.search(r"TM-score=\s*([\d.]+)\s*\(if normalized by length of Chain_1\)",
                  text)
    if m:
        result["tm_score_q"] = float(m.group(1))  # normalized by query

    m = re.search(r"TM-score=\s*([\d.]+)\s*\(if normalized by length of Chain_2\)",
                  text)
    if m:
        result["tm_score_r"] = float(m.group(1))  # normalized by reference

    return result


def parse_tmalign_alignment(text):
    """Extract the residue-level alignment mapping from TM-align stdout.

    The alignment block is at the end of stdout::

        (":" denotes aligned residue pairs of d < 5.0 A, ...)
        AG
        ::
        AG

    Returns a list of dicts::
        [
            {"query_idx": int|None, "ref_idx": int|None,
             "query_aa": str, "ref_aa": str, "distance_flag": ":" | "."},
            ...
        ]
    Where ``None`` for idx means a gap in that sequence.
    """
    lines = text.splitlines()

    # Find the line "(:" or containing the symbol legend
    legend_idx = None
    for i, line in enumerate(lines):
        if '\":\" denotes aligned' in line or \
           '\".\" denotes other' in line:
            legend_idx = i
            break

    if legend_idx is None or legend_idx + 3 > len(lines):
        return []

    # Next 3 lines (after optional blank line) are: query_seq, symbols, ref_seq
    offset = 0
    if not lines[legend_idx + 1].strip():
        offset = 1

    if legend_idx + 1 + offset + 2 >= len(lines):
        return []

    q_seq = lines[legend_idx + 1 + offset].strip()
    syms = lines[legend_idx + 2 + offset].strip()
    r_seq = lines[legend_idx + 3 + offset].strip()

    if not q_seq or not syms or not r_seq:
        return []
    if len(q_seq) != len(syms) or len(syms) != len(r_seq):
        return []

    pairs = []
    q_idx = 0
    r_idx = 0
    for qaa, sym, raa in zip(q_seq, syms, r_seq):
        if qaa != "-":
            q_idx += 1
        if raa != "-":
            r_idx += 1
        pairs.append({
            "query_idx": q_idx if qaa != "-" else None,
            "ref_idx": r_idx if raa != "-" else None,
            "query_aa": qaa,
            "ref_aa": raa,
            "distance_flag": sym,
        })

    return pairs


# ---------------------------------------------------------------------------
# Core function
# ---------------------------------------------------------------------------

def align_structures(query_path, ref_path, output_dir, skip_existing, config):
    """Align two protein structures with TM-align.

    Parameters
    ----------
    query_path  : str – path to query structure (PDB or CIF)
    ref_path    : str – path to reference structure (PDB or CIF)
    output_dir  : str – output directory
    skip_existing : bool – skip if report already exists
    config      : dict – tool config

    Returns
    -------
    dict with keys:
        report_json, aligned_pdb, scores, residue_pairs, ...
    """
    os.makedirs(output_dir, exist_ok=True)

    query_name = os.path.splitext(os.path.basename(query_path))[0]
    ref_name = os.path.splitext(os.path.basename(ref_path))[0]
    report_json = os.path.join(output_dir,
                               f"{query_name}_vs_{ref_name}.align.json")
    aligned_pdb = os.path.join(output_dir,
                               f"{query_name}_aligned.pdb")

    if skip_existing and os.path.exists(report_json):
        print(f"[SKIP] {report_json} exists, skipping alignment")
        with open(report_json) as fh:
            return json.load(fh)

    # -- convert CIF inputs to PDB if needed --
    tmp_dir = f"/tmp/align_{uuid.uuid4().hex}"
    os.makedirs(tmp_dir, exist_ok=True)

    query_pdb = query_path
    ref_pdb = ref_path

    if is_cif(query_path):
        query_pdb = os.path.join(tmp_dir, f"{query_name}_query.pdb")
        print(f"[INFO] Converting query CIF → PDB: {query_pdb}")
        cif_to_pdb(query_path, query_pdb, config)

    if is_cif(ref_path):
        ref_pdb = os.path.join(tmp_dir, f"{ref_name}_ref.pdb")
        print(f"[INFO] Converting reference CIF → PDB: {ref_pdb}")
        cif_to_pdb(ref_path, ref_pdb, config)

    # -- run TM-align --
    # Use a tmp prefix for TMalign output files (it appends suffixes)
    tm_prefix = os.path.join(tmp_dir, f"tm_{uuid.uuid4().hex}")
    cmd = (f"{config['alignment']['tmalign_bin']} {query_pdb} {ref_pdb} "
           f"-o {tm_prefix}")
    print(f"[INFO] Running TM-align: {query_name} vs {ref_name}")
    ret = run_command(cmd, "TM-align")

    # -- parse scores from stdout --
    scores = parse_tmalign_stdout(ret.stdout)

    # -- parse residue-level alignment --
    residue_pairs = parse_tmalign_alignment(ret.stdout)

    # -- collect the aligned query structure --
    # TMalign -o generates {prefix}_all_atm.pdb with the rotated/translated query
    tm_aligned = f"{tm_prefix}_all_atm.pdb"
    if os.path.exists(tm_aligned):
        shutil.copy2(tm_aligned, aligned_pdb)
        print(f"[OUT] Aligned query PDB → {aligned_pdb}")
    else:
        # Try alternative naming
        tm_aligned_alt = f"{tm_prefix}.pdb"
        if os.path.exists(tm_aligned_alt):
            shutil.copy2(tm_aligned_alt, aligned_pdb)
            print(f"[OUT] Aligned query PDB → {aligned_pdb}")
        else:
            print(f"[WARN] No aligned PDB generated by TM-align")
            aligned_pdb = None

    # -- also copy the aligned pair structure --
    aligned_pair = os.path.join(output_dir,
                                f"{query_name}_vs_{ref_name}_aligned_pair.pdb")
    for candidate in (f"{tm_prefix}_all.pdb",
                      f"{tm_prefix}_all_atm.pdb",
                      f"{tm_prefix}.pdb"):
        if os.path.exists(candidate):
            shutil.copy2(candidate, aligned_pair)
            print(f"[OUT] Aligned pair PDB → {aligned_pair}")
            break

    # -- summary dict --
    summary = {
        "query": query_name,
        "query_path": query_path,
        "reference": ref_name,
        "reference_path": ref_path,
        "scores": scores,
        "aligned_length": scores.get("aligned_length", 0),
        "rmsd": scores.get("rmsd"),
        "seq_id": scores.get("seq_id"),
        "tm_score_query": scores.get("tm_score_q"),
        "tm_score_reference": scores.get("tm_score_r"),
        "num_aligned_residues": len(residue_pairs),
        "aligned_pdb": aligned_pdb,
        "aligned_pair_pdb": aligned_pair,
        "residue_pairs": residue_pairs,
    }

    with open(report_json, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[OUT] JSON report → {report_json}")

    # -- clean up tmp --
    shutil.rmtree(tmp_dir, ignore_errors=True)

    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Align two protein structures with TM-align",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--query", "-q", required=True,
                        help="Query structure file (.pdb or .cif)")
    parser.add_argument("--reference", "-r", required=True,
                        help="Reference structure file (.pdb or .cif)")
    parser.add_argument("--output", "-o", required=True,
                        help="Output directory")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip if output report already exists")
    parser.add_argument("--config",
                        help="Path to tool config JSON (default: "
                             "data/config.json in repo root)")
    return parser.parse_args(argv)


def main():
    args = parse_args()
    config = load_config(args.config)

    if not os.path.exists(args.query):
        print(f"[ERROR] Query file not found: {args.query}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.reference):
        print(f"[ERROR] Reference file not found: {args.reference}",
              file=sys.stderr)
        sys.exit(1)

    result = align_structures(
        query_path=args.query,
        ref_path=args.reference,
        output_dir=args.output,
        skip_existing=args.skip_existing,
        config=config,
    )

    print(f"\n[DONE] Alignment complete.")
    print(f"       TM-score (ref): {result.get('tm_score_reference', 'N/A')}")
    print(f"       RMSD:           {result.get('rmsd', 'N/A')} Å")
    print(f"       Seq ID:         {result.get('seq_id', 'N/A')}")
    print(f"       Aligned length: {result.get('aligned_length', 'N/A')}")
    print(f"       JSON report:    {result.get('report_json', result.get('aligned_pdb', 'N/A'))}")


if __name__ == "__main__":
    main()
