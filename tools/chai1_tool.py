#!/usr/bin/env python3
"""
Chai-1 Structure Prediction Tool
================================

Run Chai-1 structure predictions for single or multiple sequences in a FASTA file.
Submit jobs via SLURM, wait for completion, and extract scores from predictions.

Usage:
  python chai1_tool.py --fasta sequences.fasta --output ./results \\
    --partitions "4090,3090" --max-jobs 10 --ncpus 4

Features:
  - Process single or multiple sequences from a FASTA file
  - Submit jobs to multiple SLURM partitions (evenly distribute)
  - Limit concurrent jobs with --max-jobs
  - Wait for all jobs to complete
  - Extract scores from all 5 models (pLDDT, pTM, iPTM)
  - Generate comprehensive CSV summary report
"""

import argparse
import os
import sys
import re
import subprocess as sp
import shutil
import textwrap
import csv
import time
import math
from collections import OrderedDict, defaultdict
from pathlib import Path

from _config import load_config


# ---------------------------------------------------------------------------
# Tool paths loaded from config (see data/config.json)
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# FASTA reader
# ---------------------------------------------------------------------------
def read_fasta(fpath):
    """Return list of (header, sequence) from a FASTA file."""
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


def sanitize_name(name):
    """Replace characters unsafe for filenames."""
    # Remove '>' prefix and split to get identifier
    clean = name.lstrip(">").split()[0]
    # Replace unsafe characters
    clean = re.sub(r'[/\\:<>"|?* ]', "_", clean)
    # Limit length
    return clean[:50]


# ---------------------------------------------------------------------------
# Chai-1 submission
# ---------------------------------------------------------------------------
def submit_chai1_slurm(fasta_path, output_dir, msa_path=None,
                       slurm_partition="4090", ncpus=4,
                       chai1_run=None, slurm_submit=None):
    """Submit a single Chai-1 job via SLURM.

    Returns job ID if successful, None otherwise.
    """
    # Build the command string for SLURM - run directly without cd
    chai_cmd = f"{chai1_run} {fasta_path} {output_dir}"
    if msa_path:
        chai_cmd += f" {msa_path}"

    slurm_cmd = f"{slurm_submit} \"{chai_cmd}\" {ncpus} {slurm_partition}"
    print(f"  Submitting SLURM job: {slurm_cmd}")

    result = sp.run(slurm_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    out = (result.stdout or b"").decode() + (result.stderr or b"").decode()

    # Try to extract job ID from SLURM output
    # Typical output: "Submitted batch job 123456"
    match = re.search(r"Submitted batch job (\d+)", out)
    if match:
        job_id = match.group(1)
        print(f"  Job submitted with ID: {job_id}")
        return job_id
    else:
        print(f"  WARNING: Could not extract job ID from output: {out[:200]}")
        return None


def run_chai1_local(fasta_path, output_dir, msa_path=None, chai1_run=None):
    """Run Chai-1 directly (without SLURM)."""
    cmd = [chai1_run, fasta_path, output_dir]
    if msa_path:
        cmd.append(msa_path)

    print(f"  Running: {' '.join(cmd)}")
    result = sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE)

    if result.returncode != 0:
        print(f"  WARNING: Chai-1 failed: {result.stderr.decode()[:200]}", file=sys.stderr)

    return (result.stdout or b"").decode() or (result.stderr or b"").decode()


# ---------------------------------------------------------------------------
# Job monitoring
# ---------------------------------------------------------------------------
def check_job_status(job_id):
    """Check if a SLURM job is still running.

    Returns True if job is still running/pending, False if completed/failed.
    """
    debug = False  # Set to True for debugging

    # First, try squeue (most reliable for active/pending jobs)
    cmd = f"squeue -j {job_id} --format=%T --noheader"
    result = sp.run(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)

    if result.returncode == 0:
        output = result.stdout.decode().strip()
        if debug:
            print(f"  DEBUG squeue for {job_id}: '{output}'")

        if output:
            # squeue returns states - check both full and abbreviated
            running_states = {"RUNNING", "R", "PENDING", "PD", "CONFIGURING", "CF",
                             "COMPLETING", "CG", "SUSPENDED", "S", "STOPPED", "ST"}
            running_states_lower = {s.lower() for s in running_states}
            all_running_states = running_states.union(running_states_lower)

            # Clean state (remove trailing + etc.)
            state_clean = output.rstrip('+')

            if debug:
                print(f"  DEBUG checking state: '{output}', clean: '{state_clean}'")

            if output in all_running_states or state_clean in all_running_states:
                if debug:
                    print(f"  DEBUG state in running_states -> returning True")
                return True
            else:
                # squeue returned something but not a running state
                if debug:
                    print(f"  DEBUG squeue returned non-running state: '{output}'")
                # This shouldn't happen for active jobs
                # Fall through to check sacct
                pass
        else:
            # squeue returned empty - job not in queue
            if debug:
                print(f"  DEBUG squeue returned empty")
            # Fall through to check sacct

    # If squeue doesn't show the job, check sacct (for completed/failed jobs)
    cmd = f"sacct -j {job_id} --format=State --noheader"
    result = sp.run(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)

    if result.returncode == 0:
        output = result.stdout.decode().strip()
        if debug:
            print(f"  DEBUG sacct for {job_id}: '{output}'")

        if output:
            # Check for running/pending states
            running_states = {
                "RUNNING", "R", "PENDING", "PD", "CONFIGURING", "CF",
                "COMPLETING", "CG", "SUSPENDED", "S", "STOPPED", "ST"
            }
            running_states_lower = {s.lower() for s in running_states}
            all_running_states = running_states.union(running_states_lower)

            # Check each line from sacct
            for line in output.split('\n'):
                line = line.strip()
                if not line:
                    continue

                # Get first word (state)
                parts = line.split()
                state = parts[0] if parts else ""

                # Clean state
                state_clean = state.rstrip('+')

                if debug:
                    print(f"  DEBUG checking sacct state: '{state}', clean: '{state_clean}'")

                if state in all_running_states or state_clean in all_running_states:
                    if debug:
                        print(f"  DEBUG sacct state in running_states -> returning True")
                    return True

            # If we get here, sacct didn't show any running states
            # Job is probably completed/failed
            if debug:
                print(f"  DEBUG sacct shows no running states -> returning False")
            return False
        else:
            # sacct returned empty
            if debug:
                print(f"  DEBUG sacct returned empty")
            # Job might be too new or doesn't exist
            # Be conservative: assume still running
            return True
    else:
        # sacct failed
        if debug:
            print(f"  DEBUG sacct failed")
        # Be conservative: assume still running
        return True


def wait_for_jobs(job_ids, poll_interval=30):
    """Wait for all SLURM jobs to complete."""
    if not job_ids:
        return

    print(f"\nWaiting for {len(job_ids)} jobs to complete...")
    print("Press Ctrl+C to stop waiting and proceed with completed jobs")

    remaining = set(job_ids)
    try:
        while remaining:
            completed = set()
            for job_id in remaining:
                if not check_job_status(job_id):
                    print(f"  Job {job_id} completed")
                    completed.add(job_id)

            remaining -= completed

            if remaining:
                print(f"  Still waiting for {len(remaining)} jobs...")
                time.sleep(poll_interval)
    except KeyboardInterrupt:
        print("\nInterrupted - proceeding with completed jobs")

    return list(remaining)  # Return jobs that are still running


# ---------------------------------------------------------------------------
# Score extraction
# ---------------------------------------------------------------------------
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


def load_chai_scores(npz_path):
    """Load Chai-1 NPZ scores. Returns dict with plddt, ptm, iptm, etc."""
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


def extract_all_scores(seq_name, output_dir):
    """Extract scores from all 5 models for a sequence and return only the best model.

    Returns dict with: best_model, best_plddt, best_ptm, best_iptm
    """
    scores = {}
    model_scores = []

    # Collect scores for all 5 models
    for model_idx in range(5):
        model_score = {"model_idx": model_idx}

        # Get pLDDT from CIF
        cif_path = os.path.join(output_dir, f"pred.model_idx_{model_idx}.cif")
        plddt = extract_chai_confidence(cif_path)
        if plddt is not None:
            model_score["plddt"] = plddt

        # Get pTM and iPTM from NPZ
        npz_path = os.path.join(output_dir, f"scores.model_idx_{model_idx}.npz")
        if os.path.exists(npz_path):
            npz_scores = load_chai_scores(npz_path)
            if npz_scores:
                model_score["ptm"] = npz_scores.get("ptm")
                model_score["iptm"] = npz_scores.get("iptm")

        # Only add model if we have at least pLDDT
        if "plddt" in model_score:
            model_scores.append(model_score)

    if not model_scores:
        return {}

    # Find model with highest pLDDT
    best_model = max(model_scores, key=lambda x: x["plddt"])

    # Return only best model scores
    scores["best_model"] = best_model["model_idx"]
    scores["best_plddt"] = best_model["plddt"]

    if "ptm" in best_model and best_model["ptm"] is not None:
        scores["best_ptm"] = best_model["ptm"]

    if "iptm" in best_model and best_model["iptm"] is not None:
        scores["best_iptm"] = best_model["iptm"]

    return scores


# ---------------------------------------------------------------------------
# Job distribution across partitions
# ---------------------------------------------------------------------------
def distribute_jobs_across_partitions(jobs, partitions, max_jobs=None):
    """Distribute jobs evenly across available partitions.

    Args:
        jobs: List of job tuples (seq_name, fasta_path, output_dir)
        partitions: List of partition names
        max_jobs: Maximum number of jobs to submit at once

    Returns:
        List of (partition, job) assignments
    """
    if max_jobs and len(jobs) > max_jobs:
        print(f"Limiting to {max_jobs} jobs (from {len(jobs)})")
        jobs = jobs[:max_jobs]

    assignments = []
    for i, job in enumerate(jobs):
        partition = partitions[i % len(partitions)]
        assignments.append((partition, job))

    return assignments


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------
def write_summary_report(sequences, all_scores, output_dir):
    """Write CSV report with best model scores for each sequence."""
    report_path = os.path.join(output_dir, "chai1_scores_summary.csv")

    # Define fixed column order
    fieldnames = [
        "seq_name", "seq_header", "sequence", "length",
        "best_model", "best_plddt", "best_ptm", "best_iptm"
    ]

    with open(report_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for (header, seq), scores in zip(sequences, all_scores):
            seq_name = sanitize_name(header)
            row = {
                "seq_name": seq_name,
                "seq_header": header.lstrip(">"),
                "sequence": seq,
                "length": len(seq),
                "best_model": scores.get("best_model", ""),
                "best_plddt": f"{scores.get('best_plddt', ''):.2f}" if scores.get('best_plddt') else "",
                "best_ptm": f"{scores.get('best_ptm', ''):.4f}" if scores.get('best_ptm') else "",
                "best_iptm": f"{scores.get('best_iptm', ''):.4f}" if scores.get('best_iptm') else "",
            }
            writer.writerow(row)

    print(f"\nSummary report written to: {report_path}")
    return report_path


def print_summary_table(sequences, all_scores):
    """Print a formatted summary table to console."""
    print("\n" + "=" * 80)
    print("CHAI-1 PREDICTION SUMMARY (Best Model Only)")
    print("=" * 80)
    print(f"{'Sequence':<25} {'Length':<8} {'Model':<6} {'pLDDT':<8} {'pTM':<10} {'iPTM':<10}")
    print("-" * 80)

    for (header, seq), scores in zip(sequences, all_scores):
        seq_name = sanitize_name(header)
        short_name = seq_name[:23] + "..." if len(seq_name) > 23 else seq_name
        best_model = scores.get("best_model", "N/A")
        best_plddt = scores.get("best_plddt", "N/A")
        best_ptm = scores.get("best_ptm", "N/A")
        best_iptm = scores.get("best_iptm", "N/A")

        if isinstance(best_plddt, float):
            best_plddt = f"{best_plddt:.2f}"
        if isinstance(best_ptm, float):
            best_ptm = f"{best_ptm:.4f}"
        if isinstance(best_iptm, float):
            best_iptm = f"{best_iptm:.4f}"

        print(f"{short_name:<25} {len(seq):<8} {best_model:<6} {best_plddt:<8} {best_ptm:<10} {best_iptm:<10}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Chai-1 Structure Prediction Tool - Run predictions for single or multiple sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples:
              # Single sequence with default partition
              %(prog)s --fasta sequence.fasta --output ./results

              # Multiple sequences, distribute across partitions
              %(prog)s --fasta sequences.fasta --output ./results \\
                --partitions "4090,3090" --max-jobs 20 --ncpus 4

              # With custom MSA for each sequence
              %(prog)s --fasta sequences.fasta --output ./results \\
                --msa-dir ./msa_files --msa-suffix ".a3m"

              # Skip SLURM, run locally (for testing)
              %(prog)s --fasta sequence.fasta --output ./results --local
        """),
    )

    ap.add_argument("--fasta", required=True, help="Input FASTA file (single or multiple sequences)")
    ap.add_argument("--output", "-o", default="./chai1_results",
                   help="Output directory (default: ./chai1_results)")
    ap.add_argument("--partitions", default="4090",
                   help="Comma-separated list of SLURM partitions (default: 4090)")
    ap.add_argument("--max-jobs", type=int, default=None,
                   help="Maximum number of jobs to submit at once (default: no limit)")
    ap.add_argument("--ncpus", type=int, default=4,
                   help="Number of CPUs per job (default: 4)")
    ap.add_argument("--msa-dir",
                   help="Directory containing MSA files (named after sequence IDs)")
    ap.add_argument("--msa-suffix", default=".a3m",
                   help="Suffix for MSA files (default: .a3m)")
    ap.add_argument("--local", action="store_true",
                   help="Run Chai-1 locally instead of submitting to SLURM")
    ap.add_argument("--no-wait", action="store_true",
                   help="Don't wait for jobs to complete (submit and exit)")
    ap.add_argument("--poll-interval", type=int, default=30,
                   help="Seconds between job status checks (default: 30)")
    ap.add_argument("--skip-existing", action="store_true",
                   help="Skip sequences that already have prediction results")
    ap.add_argument("--config", default=None,
                   help="Path to config JSON file (default: <repo_root>/data/config.json)")

    args = ap.parse_args()

    # Load tool dependency config
    cfg = load_config(args.config)
    chai1_dir = cfg["chai1"]["chai1_dir"]
    chai1_run = os.path.join(chai1_dir, "run.sh")
    slurm_submit = cfg["slurm"]["submit_script"]

    # -----------------------------------------------------------------------
    # Setup
    # -----------------------------------------------------------------------
    print("=" * 60)
    print("CHAI-1 STRUCTURE PREDICTION TOOL")
    print("=" * 60)

    root = os.path.abspath(args.output)
    os.makedirs(root, exist_ok=True)

    partitions = [p.strip() for p in args.partitions.split(",") if p.strip()]
    print(f"Output directory: {root}")
    print(f"SLURM partitions: {partitions}")
    print(f"CPUs per job: {args.ncpus}")
    print(f"Local mode: {args.local}")

    # -----------------------------------------------------------------------
    # Read sequences
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 1: Reading sequences")
    print("=" * 60)

    sequences = read_fasta(args.fasta)
    if not sequences:
        print(f"ERROR: No sequences found in {args.fasta}")
        sys.exit(1)

    print(f"Found {len(sequences)} sequence(s)")
    for i, (header, seq) in enumerate(sequences[:5]):  # Show first 5
        print(f"  {i+1}. {header[:50]}{'...' if len(header) > 50 else ''} ({len(seq)} aa)")
    if len(sequences) > 5:
        print(f"  ... and {len(sequences) - 5} more")

    # -----------------------------------------------------------------------
    # Prepare jobs
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 2: Preparing jobs")
    print("=" * 60)

    # Create fasta directory for input files
    fasta_dir = os.path.join(root, "fastas")
    os.makedirs(fasta_dir, exist_ok=True)

    jobs = []
    existing_results = []  # Track sequences that already have results

    for header, seq in sequences:
        seq_name = sanitize_name(header)

        # Create FASTA file in fasta directory
        fasta_path = os.path.join(fasta_dir, f"{seq_name}.fasta")

        # Output directory for Chai-1 (will be created by Chai-1)
        output_dir = os.path.join(root, seq_name)

        # Check if results already exist
        has_existing_results = False
        if args.skip_existing:
            # Check if all 5 models exist
            all_exist = True
            for model_idx in range(5):
                cif_path = os.path.join(output_dir, f"pred.model_idx_{model_idx}.cif")
                if not os.path.exists(cif_path):
                    all_exist = False
                    break

            if all_exist:
                print(f"  {seq_name}: prediction results already exist (will skip submission)")
                has_existing_results = True
                existing_results.append(seq_name)

        # Write FASTA file (always write it, even for existing results)
        with open(fasta_path, "w") as fh:
            fh.write(f">protein|{header.lstrip('>')}\n{seq}\n")

        # Find MSA file if provided
        msa_path = None
        if args.msa_dir:
            # Try different naming patterns
            msa_patterns = [
                os.path.join(args.msa_dir, f"{seq_name}{args.msa_suffix}"),
                os.path.join(args.msa_dir, f"{seq_name.replace('_', '')}{args.msa_suffix}"),
                os.path.join(args.msa_dir, f"{header.lstrip('>').split()[0]}{args.msa_suffix}"),
            ]
            for pattern in msa_patterns:
                if os.path.exists(pattern):
                    msa_path = pattern
                    break

        # Add to jobs list with a flag indicating if results already exist
        jobs.append((seq_name, fasta_path, output_dir, msa_path, has_existing_results))

    if not jobs:
        print("No jobs to process (all skipped or no sequences)")
        sys.exit(0)

    print(f"Prepared {len(jobs)} job(s) for processing")

    # -----------------------------------------------------------------------
    # Submit jobs
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 3: Submitting jobs")
    print("=" * 60)

    if args.local:
        # Run locally
        jobs_to_submit = [(seq_name, fasta_path, output_dir, msa_path)
                         for seq_name, fasta_path, output_dir, msa_path, has_existing in jobs
                         if not has_existing]

        if jobs_to_submit:
            for seq_name, fasta_path, output_dir, msa_path in jobs_to_submit:
                print(f"\nRunning Chai-1 for {seq_name}...")
                run_chai1_local(fasta_path, output_dir, msa_path, chai1_run=chai1_run)
            print(f"\nSubmitted {len(jobs_to_submit)} job(s) locally")
        else:
            print("\nNo jobs to submit (all sequences already have results)")
    else:
        # Filter out jobs that already have results
        jobs_to_submit = [(seq_name, fasta_path, output_dir, msa_path)
                         for seq_name, fasta_path, output_dir, msa_path, has_existing in jobs
                         if not has_existing]

        if jobs_to_submit:
            # Distribute jobs across partitions
            assignments = distribute_jobs_across_partitions(jobs_to_submit, partitions, args.max_jobs)

            # Submit jobs and collect job IDs
            job_ids = []
            partition_counts = defaultdict(int)

            for partition, (seq_name, fasta_path, output_dir, msa_path) in assignments:
                print(f"\nSubmitting {seq_name} to partition {partition}...")
                job_id = submit_chai1_slurm(fasta_path, output_dir, msa_path, partition, args.ncpus,
                                            chai1_run=chai1_run, slurm_submit=slurm_submit)

                if job_id:
                    job_ids.append(job_id)
                    partition_counts[partition] += 1

            print(f"\nSubmitted {len(job_ids)} job(s) total:")
            for partition, count in partition_counts.items():
                print(f"  {partition}: {count} job(s)")
        else:
            print("\nNo jobs to submit (all sequences already have results)")

        # Wait for jobs to complete
        if not args.no_wait and job_ids:
            still_running = wait_for_jobs(job_ids, args.poll_interval)
            if still_running:
                print(f"\nWarning: {len(still_running)} jobs still running")
        elif args.no_wait:
            print("\nNot waiting for jobs to complete (--no-wait flag)")

    # -----------------------------------------------------------------------
    # Extract scores
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 4: Extracting scores")
    print("=" * 60)

    all_scores = []
    successful = 0

    for seq_name, _, output_dir, _, _ in jobs:
        print(f"\nExtracting scores for {seq_name}...")

        # Check if prediction completed
        has_models = False
        for model_idx in range(5):
            cif_path = os.path.join(output_dir, f"pred.model_idx_{model_idx}.cif")
            if os.path.exists(cif_path):
                has_models = True
                break

        if not has_models:
            print(f"  WARNING: No prediction results found for {seq_name}")
            all_scores.append({})
            continue

        # Extract scores
        scores = extract_all_scores(seq_name, output_dir)
        all_scores.append(scores)

        if scores:
            successful += 1
            # Print key scores
            best_model = scores.get("best_model", "N/A")
            best_plddt = scores.get("best_plddt", "N/A")
            best_ptm = scores.get("best_ptm", "N/A")

            if isinstance(best_plddt, float):
                print(f"  Best model: {best_model}, pLDDT: {best_plddt:.2f}, pTM: {best_ptm:.4f}")
            else:
                print(f"  Scores extracted but some values missing")
        else:
            print(f"  No scores could be extracted")

    print(f"\nSuccessfully extracted scores for {successful}/{len(jobs)} sequence(s)")

    # -----------------------------------------------------------------------
    # Generate reports
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("STEP 5: Generating reports")
    print("=" * 60)

    # Write detailed CSV report
    # all_scores has one entry per sequence (same as sequences)
    report_path = write_summary_report(sequences, all_scores, root)

    # Print summary table
    print_summary_table(sequences, all_scores)

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)
    print(f"Output directory: {root}")
    print(f"  FASTA files:       {fasta_dir}")
    print(f"  Prediction results: {root}/<sequence_name>/")
    print(f"Total sequences:     {len(sequences)}")
    if args.skip_existing and existing_results:
        print(f"  - With existing results: {len(existing_results)} (skipped submission)")
        print(f"  - New predictions needed: {len(jobs) - len(existing_results)}")
    print(f"Scores extracted:    {successful}/{len(jobs)}")
    print(f"Detailed report:     {report_path}")
    print(f"\nTo view top predictions by pLDDT:")
    print(f"  sort -t, -k1,1 -k7,7nr {report_path} | head -20")


if __name__ == "__main__":
    main()