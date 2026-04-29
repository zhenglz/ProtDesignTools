#!/usr/bin/env python3
"""
GOAP Energy Scoring Tool
=========================

Calculate protein statistical energy using the Fast_GOAP program.
Accepts a single PDB/CIF file or a directory of PDB/CIF files.
CIF files are automatically converted to PDB using cif2pdb.py.

Fast_GOAP outputs three tab-separated values per model:
  GOAP-AG  (atom-general pairwise potential)
  GOAP-AS  (atom-specific potential)
  GOAP-Total (total statistical potential energy, lower is better)

Outputs:
  1. CSV file with all GOAP scores sorted by total energy
  2. Plain-text summary

Tools and paths are defined in data/config.json under "goap".

Usage:
  python tools/goap_tool.py -i structure.pdb -o ./goap_results
  python tools/goap_tool.py -i structure.cif -o ./goap_results
  python tools/goap_tool.py -d ./pdb_dir -o ./goap_results
  python tools/goap_tool.py -d ./cif_dir -o ./goap_results --suffix .cif
"""

import argparse
import os
import shutil
import subprocess as sp
import sys
import tempfile
import uuid

import pandas as pd

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


def is_cif(fpath):
    """Return True if the file is an mmCIF format."""
    ext = os.path.splitext(fpath)[1].lstrip(".").lower()
    return ext in ("cif", "mmcif")


def cif_to_pdb(cif_path, pdb_path, config):
    """Convert a CIF file to PDB using cif2pdb.py."""
    python_exe = config["alignment"]["python_exe"]
    script = config["alignment"]["cif2pdb_script"]
    cmd = f"{python_exe} {script} {cif_path} {pdb_path} --quiet"
    run_command(cmd, f"cif_to_pdb: {os.path.basename(cif_path)}")
    if not os.path.exists(pdb_path) or os.path.getsize(pdb_path) == 0:
        print(f"[ERROR] CIF->PDB conversion produced empty file: {pdb_path}",
              file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# GOAP output parser
# ---------------------------------------------------------------------------

def parse_goap_stdout(text):
    """Parse Fast_GOAP tab-separated stdout.

    Fast_GOAP outputs one line per model::

        model_name\tgoap_ag\tgoap_as\tgoap_total

    Returns list of dicts with keys: model, goap_ag, goap_as, goap_total.
    """
    results = []
    for line in text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 4:
            try:
                results.append({
                    "model": parts[0],
                    "goap_ag": float(parts[1]),
                    "goap_as": float(parts[2]),
                    "goap_total": float(parts[3]),
                })
            except ValueError:
                continue
    return results


# ---------------------------------------------------------------------------
# Core function
# ---------------------------------------------------------------------------

def run_goap(input_path, output_dir, suffix, skip_existing, config):
    """Run Fast_GOAP on a single PDB file or a directory of PDB files.

    Parameters
    ----------
    input_path  : str – path to a single PDB/CIF file, or a directory
    output_dir  : str – output directory
    suffix      : str – file suffix for directory mode ('.pdb' or '.cif')
    skip_existing : bool – skip if CSV output already exists
    config      : dict – tool config

    Returns
    -------
    dict with keys: csv_path, results (list of dicts), num_models
    """
    os.makedirs(output_dir, exist_ok=True)
    goap_bin = config["goap"]["goap_bin"]

    csv_path = os.path.join(output_dir, "goap_scores.csv")
    summary_path = os.path.join(output_dir, "goap_summary.txt")

    if skip_existing and os.path.exists(csv_path):
        print(f"[SKIP] {csv_path} exists, skipping GOAP calculation")
        df = pd.read_csv(csv_path)
        return {
            "csv_path": csv_path,
            "summary_path": summary_path,
            "results": df.to_dict("records"),
            "num_models": len(df),
        }

    is_single_file = os.path.isfile(input_path)
    tmp_dir = f"/tmp/goap_{uuid.uuid4().hex}"

    if is_single_file:
        # -- Single file mode --
        if is_cif(input_path):
            os.makedirs(tmp_dir, exist_ok=True)
            pdb_name = os.path.splitext(os.path.basename(input_path))[0] + ".pdb"
            work_pdb = os.path.join(tmp_dir, pdb_name)
            cif_to_pdb(input_path, work_pdb, config)
        else:
            work_pdb = input_path

        print(f"[INFO] Running Fast_GOAP on: {os.path.basename(work_pdb)}")
        ret = sp.run([goap_bin, "-i", work_pdb],
                     stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
        if ret.returncode != 0:
            print(f"[ERROR] Fast_GOAP failed (exit {ret.returncode})",
                  file=sys.stderr)
            if ret.stderr:
                print(ret.stderr, file=sys.stderr)
            sys.exit(1)

        results = parse_goap_stdout(ret.stdout)

    else:
        # -- Directory mode --
        is_cif_suffix = suffix.lstrip(".").lower() in ("cif", "mmcif")

        if is_cif_suffix:
            # Convert all CIF files to PDB in a temp directory
            os.makedirs(tmp_dir, exist_ok=True)
            cif_files = sorted([
                f for f in os.listdir(input_path)
                if f.endswith(suffix)
            ])
            if not cif_files:
                print(f"[ERROR] No '{suffix}' files found in {input_path}",
                      file=sys.stderr)
                sys.exit(1)
            print(f"[INFO] Converting {len(cif_files)} CIF file(s) to PDB...")
            for cif_file in cif_files:
                cif_fpath = os.path.join(input_path, cif_file)
                pdb_name = os.path.splitext(cif_file)[0] + ".pdb"
                pdb_fpath = os.path.join(tmp_dir, pdb_name)
                cif_to_pdb(cif_fpath, pdb_fpath, config)
            effective_dir = tmp_dir
            effective_suffix = ".pdb"
        else:
            effective_dir = input_path
            effective_suffix = suffix

        # Build list file (base names without suffix for Fast_GOAP list mode)
        pdb_files = sorted([
            f for f in os.listdir(effective_dir)
            if f.endswith(effective_suffix)
        ])
        if not pdb_files:
            print(f"[ERROR] No '{effective_suffix}' files found in "
                  f"{effective_dir}", file=sys.stderr)
            sys.exit(1)

        # Fast_GOAP list mode: entries are basenames without suffix
        list_fd, list_path = tempfile.mkstemp(suffix=".list", text=True)
        with os.fdopen(list_fd, "w") as fh:
            for f in pdb_files:
                base = f[:-len(effective_suffix)]  # strip suffix
                fh.write(base + "\n")

        print(f"[INFO] Running Fast_GOAP on {len(pdb_files)} file(s) "
              f"in {effective_dir}")
        ret = sp.run([goap_bin, "-L", list_path,
                      "-r", effective_dir, "-s", effective_suffix],
                     stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
        os.unlink(list_path)

        if ret.returncode != 0:
            print(f"[ERROR] Fast_GOAP failed (exit {ret.returncode})",
                  file=sys.stderr)
            if ret.stderr:
                print(ret.stderr, file=sys.stderr)
            sys.exit(1)

        results = parse_goap_stdout(ret.stdout)

    # -- clean up temp dir --
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir, ignore_errors=True)

    if not results:
        print("[WARN] No GOAP scores parsed from Fast_GOAP output")
        return {
            "csv_path": None,
            "summary_path": None,
            "results": [],
            "num_models": 0,
        }

    # -- sort by total energy (lower is better) --
    results.sort(key=lambda x: x["goap_total"])

    # -- write CSV --
    df = pd.DataFrame(results)
    df.to_csv(csv_path, index=False, float_format="%.6f")
    print(f"[OUT] GOAP scores CSV -> {csv_path}")

    # -- write text summary --
    with open(summary_path, "w") as fh:
        fh.write("GOAP Statistical Energy Scores\n")
        fh.write("==============================\n")
        if is_single_file:
            fh.write(f"Input: {input_path}\n")
        else:
            fh.write(f"Input directory: {input_path}\n")
            fh.write(f"Suffix: {suffix}\n")
        fh.write(f"Number of models: {len(results)}\n\n")
        fh.write("Rank  Model                          "
                 "GOAP-AG         GOAP-AS         GOAP-Total\n")
        fh.write("-" * 90 + "\n")
        for i, r in enumerate(results, 1):
            fh.write(f"{i:4d}  {r['model']:<30s}  "
                     f"{r['goap_ag']:14.6f}  {r['goap_as']:14.6f}  "
                     f"{r['goap_total']:14.6f}\n")
    print(f"[OUT] GOAP summary -> {summary_path}")

    return {
        "csv_path": csv_path,
        "summary_path": summary_path,
        "results": results,
        "num_models": len(results),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Calculate protein statistical energy using Fast_GOAP",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    # Input selection (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-i", "--input",
        help="Input PDB or CIF file (single-file mode)")
    input_group.add_argument(
        "-d", "--dir",
        help="Input directory containing PDB or CIF files (batch mode)")

    parser.add_argument("-o", "--output", required=True,
                        help="Output directory for GOAP results")
    parser.add_argument("-s", "--suffix", default=".pdb",
                        help="File suffix for directory mode "
                             "(default: .pdb; use .cif for CIF files)")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip calculation if output CSV already exists")
    parser.add_argument("--config",
                        help="Path to tool config JSON (default: "
                             "data/config.json in repo root)")
    return parser.parse_args(argv)


def main():
    args = parse_args()
    config = load_config(args.config)

    # Validate GOAP configuration
    if "goap" not in config or "goap_bin" not in config.get("goap", {}):
        print("[ERROR] 'goap.goap_bin' not found in config.json",
              file=sys.stderr)
        sys.exit(1)

    goap_bin = config["goap"]["goap_bin"]
    if not os.path.exists(goap_bin):
        print(f"[ERROR] Fast_GOAP binary not found: {goap_bin}",
              file=sys.stderr)
        sys.exit(1)

    # Determine input path
    input_path = args.input if args.input else args.dir
    if not os.path.exists(input_path):
        print(f"[ERROR] Input path not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    result = run_goap(
        input_path=input_path,
        output_dir=args.output,
        suffix=args.suffix,
        skip_existing=args.skip_existing,
        config=config,
    )

    if result["num_models"] == 0:
        print("\n[DONE] No GOAP scores computed.")
        return

    # Print top-10 summary to stdout
    print(f"\n[DONE] GOAP scoring complete.")
    print(f"       Models scored: {result['num_models']}")
    print(f"       CSV:           {result['csv_path']}")
    print(f"       Summary:       {result['summary_path']}")
    print(f"\nTop 10 (lower total energy = better):")
    print(f"{'Rank':<6s}{'Model':<32s}{'GOAP-AG':>14s}{'GOAP-AS':>14s}"
          f"{'GOAP-Total':>14s}")
    print("-" * 80)
    for i, r in enumerate(result["results"][:10], 1):
        print(f"{i:<6d}{r['model']:<32s}{r['goap_ag']:14.6f}"
              f"{r['goap_as']:14.6f}{r['goap_total']:14.6f}")


if __name__ == "__main__":
    main()
