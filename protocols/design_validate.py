#!/usr/bin/env python3
"""
Protocol: RFDiffusion3 Backbone Design → Chai-1 Structure Validation.

Pipeline:
  1. RFDiffusion3: generate N backbone designs for specified region(s)
  2. Select top M by composite score
  3. Chai-1: predict structures for top M designs (SLURM)
  4. Rank by pLDDT * iPTM, select top K

Usage:
  python protocols/design_validate.py \\
    --pdb input.pdb --design-regions "A90-95:5-5" \\
    --num-designs 50 --top-m 10 --top-k 3 --output ./results
"""

import argparse, csv, os, sys, subprocess as sp, time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / 'tools'))

from _config import load_config
from chai1_tool import submit_chai1_slurm, check_job_status, extract_all_scores

RFD3_SCRIPT = str(REPO_ROOT / 'tools' / 'rfdiffusion3_tool.py')


# ---------------------------------------------------------------------------
# Composite rank (same formula as rfdiffusion3_tool.py — lower = better)
# ---------------------------------------------------------------------------
def composite_rank(d):
    return (
        int(d.get('n_chainbreaks', 0)) * 100 +
        int(d.get('n_clashing', 0)) * 10 +
        abs(float(d.get('radius_of_gyration', 0)) - 12) * 0.1 +
        float(d.get('alanine_content', 0)) * 10 +
        abs(float(d.get('glycine_content', 0)) - 0.05) * 10 +
        float(d.get('max_ca_deviation', 0))
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description='RFDiffusion3 Design → Chai-1 Validation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument('--pdb', required=True, help='Input PDB file')
    p.add_argument('--design-regions', required=True,
                   help='Design regions: <chain><start>-<end>:<len_min>-<len_max>[;...]')
    p.add_argument('--num-designs', type=int, default=50,
                   help='RFD3 designs to generate (default: 50)')
    p.add_argument('--top-m', type=int, default=10,
                   help='Top designs to validate with Chai-1 (default: 10)')
    p.add_argument('--top-k', type=int, default=3,
                   help='Final top designs to report (default: 3)')
    p.add_argument('--output', '-o', default='./protocol_out')
    p.add_argument('--slurm-partition', default='4090')
    p.add_argument('--ncpus', type=int, default=4)
    p.add_argument('--max-jobs', type=int, default=4)
    p.add_argument('--local', action='store_true',
                   help='Run locally instead of SLURM')
    p.add_argument('--config', help='Path to config JSON')
    return p.parse_args()


# ---------------------------------------------------------------------------
# Phase helpers
# ---------------------------------------------------------------------------

def run_rfd3(args, rfd3_out, config_path):
    """Run RFDiffusion3 as subprocess and return output directory."""
    cmd = [
        sys.executable, RFD3_SCRIPT,
        '--pdb', args.pdb,
        '--design-regions', args.design_regions,
        '--num-designs', str(args.num_designs),
        '--output', rfd3_out,
        '--top-n', '0',
        '--max-jobs', str(args.max_jobs),
    ]
    if config_path:
        cmd += ['--config', config_path]
    if args.local:
        cmd.append('--local')
    else:
        cmd += ['--slurm-partition', args.slurm_partition,
                '--ncpus', str(args.ncpus)]

    print(f'[Phase 1] Running RFDiffusion3 ({args.num_designs} designs)...')
    sp.run(cmd, check=True)
    return rfd3_out


def load_rfd3_scores(rfd3_out):
    """Read rfd3_scores.csv, return list of dicts with composite score."""
    csv_path = os.path.join(rfd3_out, 'rfd3_scores.csv')
    if not os.path.exists(csv_path):
        print(f'ERROR: scores CSV not found at {csv_path}')
        sys.exit(1)
    designs = []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            row['_score'] = composite_rank(row)
            designs.append(row)
    return designs


def run_chai1_batch(designs, rfd3_out, chai1_out, cfg, args):
    """Submit Chai-1 SLURM jobs for each design, wait, extract scores."""
    os.makedirs(chai1_out, exist_ok=True)
    chai1_run = os.path.join(cfg['chai1']['chai1_dir'], 'run.sh')
    slurm_submit = cfg['slurm']['submit_script']
    fasta_dir = os.path.join(rfd3_out, 'top_designs_fastas')

    print(f'[Phase 3] Submitting Chai-1 jobs for {len(designs)} designs...')

    # Submit all jobs
    jobs = {}  # job_id -> (name, out_dir)
    for d in designs:
        name = d['design_name']
        fasta = os.path.join(fasta_dir, f'{name}.fasta')
        if not os.path.exists(fasta):
            print(f'  SKIP {name}: no FASTA found')
            continue
        out_dir = os.path.join(chai1_out, name)
        jid = submit_chai1_slurm(
            fasta, out_dir,
            slurm_partition=args.slurm_partition, ncpus=args.ncpus,
            chai1_run=chai1_run, slurm_submit=slurm_submit,
        )
        if jid:
            jobs[jid] = (name, out_dir)
            print(f'  {name}: submitted (job {jid})')

    if not jobs:
        print('No jobs submitted.')
        return []

    # Wait for all jobs
    total = len(jobs)
    while jobs:
        done = []
        for jid, (name, out_dir) in jobs.items():
            st = check_job_status(jid)
            if st in ('COMPLETED', 'CD', 'COMPLETING'):
                done.append(jid)
            elif st in ('FAILED', 'F', 'CANCELLED', 'CA', 'TIMEOUT', 'TO'):
                print(f'  {name}: FAILED ({st})')
                done.append(jid)
        for jid in done:
            del jobs[jid]
        if jobs:
            print(f'  [{time.strftime("%H:%M:%S")}] {total - len(jobs)}/{total} done, '
                  f'{len(jobs)} pending')
            time.sleep(30)

    # Extract scores
    results = []
    for d in designs:
        name = d['design_name']
        out_dir = os.path.join(chai1_out, name)
        scores = extract_all_scores(name, out_dir)
        if scores:
            results.append({
                'design_name': name,
                'plddt': scores.get('best_plddt', 0),
                'ptm': scores.get('best_ptm', 0),
                'iptm': scores.get('best_iptm', 0),
                'combined': scores.get('best_plddt', 0) * scores.get('best_iptm', 0),
            })
    return results


def write_final_results(results, output_dir):
    """Save final top-K results to CSV and print summary."""
    csv_path = os.path.join(output_dir, 'final_topk.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['rank', 'design_name', 'plddt', 'ptm', 'iptm', 'combined'])
        for i, r in enumerate(results, 1):
            w.writerow([i, r['design_name'], r['plddt'], r['ptm'], r['iptm'], r['combined']])
    print(f'\nResults saved to {csv_path}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    cfg = load_config(args.config)
    os.makedirs(args.output, exist_ok=True)

    rfd3_out = os.path.join(args.output, 'rfd3_designs')

    # Phase 1: RFDiffusion3 backbone generation
    run_rfd3(args, rfd3_out, args.config)

    # Phase 2: Select top M designs by composite score
    designs = load_rfd3_scores(rfd3_out)
    designs.sort(key=lambda d: d['_score'])
    top_m = designs[:args.top_m]
    print(f'[Phase 2] Selected top {len(top_m)}/{len(designs)} designs:')
    for i, d in enumerate(top_m):
        print(f'  {i+1}. {d["design_name"]}  score={d["_score"]:.1f}')

    if not top_m:
        print('No designs selected. Exiting.')
        return

    # Phase 3: Chai-1 prediction + Phase 4: Final ranking
    results = run_chai1_batch(top_m, rfd3_out,
                              os.path.join(args.output, 'chai1_preds'),
                              cfg, args)
    results.sort(key=lambda r: r['combined'], reverse=True)
    top_k = results[:args.top_k]

    print(f'\n[Phase 4] Top {len(top_k)} final designs:')
    for i, r in enumerate(top_k):
        print(f'  {i+1}. {r["design_name"]:30s} '
              f'pLDDT={r["plddt"]:.3f}  pTM={r["ptm"]:.3f}  iPTM={r["iptm"]:.3f}  '
              f'combined={r["combined"]:.3f}')

    write_final_results(top_k, args.output)
    print('\nDone.')


if __name__ == '__main__':
    main()
