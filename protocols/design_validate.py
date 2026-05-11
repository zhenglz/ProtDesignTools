#!/usr/bin/env python3
"""
Protocol: RFDiffusion3 Backbone Design → Chai-1 Structure Validation.

Pipeline:
  1. RFDiffusion3: generate N backbone designs across multiple rounds
  2. Collect all designs, deduplicate by designed-region sequences
  3. Select top M unique designs by composite score
  4. Chai-1: predict structures for top M designs (SLURM)
  5. Rank by pLDDT * iPTM, select top K

Usage:
  python protocols/design_validate.py \\
    --pdb input.pdb --design-regions "A90-95:5-5" \\
    --num-rounds 3 --num-designs 50 --top-m 20 --top-k 5 --output ./results
"""

import argparse, csv, os, sys, subprocess as sp, time, shutil
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / 'tools'))

from _config import load_config
from chai1_tool import submit_chai1_slurm, check_job_status, extract_all_scores
from design_utils import (
    composite_rank, collect_designs_from_rounds,
    deduplicate_designs, write_unique_fastas,
)

RFD3_SCRIPT = str(REPO_ROOT / 'tools' / 'rfdiffusion3_tool.py')


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
    p.add_argument('--num-rounds', type=int, default=1,
                   help='Number of independent RFD3 rounds (default: 1)')
    p.add_argument('--num-designs', type=int, default=50,
                   help='RFD3 designs per round (default: 50)')
    p.add_argument('--top-m', type=int, default=10,
                   help='Top unique designs to validate with Chai-1 (default: 10)')
    p.add_argument('--top-k', type=int, default=3,
                   help='Final top designs to report (default: 3)')
    p.add_argument('--output', '-o', default='./protocol_out')
    p.add_argument('--slurm-partition', default='4090')
    p.add_argument('--ncpus', type=int, default=4)
    p.add_argument('--max-jobs', type=int, default=4)
    p.add_argument('--local', action='store_true',
                   help='Run locally instead of SLURM')
    p.add_argument('--config', help='Path to config JSON')
    p.add_argument('--skip-existing', action='store_true',
                   help='Skip phases where output files already exist')
    return p.parse_args()


# ---------------------------------------------------------------------------
# Phase helpers
# ---------------------------------------------------------------------------

def build_rfd3_cmd(args, rfd3_out, config_path):
    """Build the RFDiffusion3 subprocess command list."""
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
    return cmd


def run_chai1_batch(designs, fasta_dir, chai1_out, cfg, args):
    """Submit Chai-1 SLURM jobs for each design, wait, extract scores.

    Args:
        designs: list of design dicts, each with '_tag' and '_unique_fasta'
        fasta_dir: directory containing per-design FASTA files
        chai1_out: output directory for Chai-1 predictions
        cfg: config dict
        args: CLI args (for slurm_partition, ncpus, local, skip_existing)
    """
    os.makedirs(chai1_out, exist_ok=True)
    chai1_run = os.path.join(cfg['chai1']['chai1_dir'], 'run.sh')
    slurm_submit = cfg['slurm']['submit_script']

    # Separate into submit-vs-skip
    to_submit = []
    skipped = []
    for d in designs:
        tag = d.get('_tag', d['design_name'])
        out_dir = os.path.join(chai1_out, tag)
        # Check all 5 model outputs exist
        all_five = all(
            os.path.exists(os.path.join(out_dir, f'pred.model_idx_{m}.cif'))
            for m in range(5)
        )
        if args.skip_existing and all_five:
            skipped.append((tag, out_dir, d))
        else:
            to_submit.append(d)

    if skipped:
        print(f'  Skipping {len(skipped)} design(s) with existing results')

    # Submit jobs
    jobs = {}  # job_id -> (tag, out_dir)
    if to_submit:
        print(f'[{time.strftime("%H:%M:%S")}] Submitting {len(to_submit)} Chai-1 job(s)...')
        for d in to_submit:
            tag = d.get('_tag', d['design_name'])
            fasta = d.get('_unique_fasta', os.path.join(fasta_dir, f'{tag}.fasta'))
            if not os.path.exists(fasta):
                print(f'  SKIP {tag}: no FASTA found')
                continue
            out_dir = os.path.join(chai1_out, tag)
            jid = submit_chai1_slurm(
                fasta, out_dir,
                slurm_partition=args.slurm_partition, ncpus=args.ncpus,
                chai1_run=chai1_run, slurm_submit=slurm_submit,
            )
            if jid:
                jobs[jid] = (tag, out_dir)
                print(f'  {tag}: submitted (job {jid})')

    if not jobs and not skipped:
        print('No jobs submitted and no existing results.')
        return []

    # Wait for jobs that were submitted
    total = len(jobs)
    while jobs:
        done = []
        for jid, (tag, out_dir) in jobs.items():
            still_running = check_job_status(jid)  # True=running/pending, False=completed/failed
            if not still_running:
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
        tag = d.get('_tag', d['design_name'])
        out_dir = os.path.join(chai1_out, tag)
        scores = extract_all_scores(tag, out_dir)
        if scores:
            results.append({
                'design_name': tag,
                'orig_name': d['design_name'],
                'best_model': scores.get('best_model', ''),
                'plddt': scores.get('best_plddt', 0),
                'ptm': scores.get('best_ptm', 0),
                'iptm': scores.get('best_iptm', 0),
                'combined': scores.get('best_plddt', 0) * scores.get('best_iptm', 0),
            })
    return results


def _cif_to_pdb(cif_path, pdb_path):
    """Convert a CIF file to PDB format using BioPython."""
    try:
        from Bio.PDB import MMCIFParser, PDBIO
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("model", cif_path)
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_path)
        return True
    except ImportError:
        print("  BioPython not available, copying CIF as-is instead")
        shutil.copy2(cif_path, pdb_path)
        return False
    except Exception as e:
        print(f"  CIF→PDB conversion failed: {e}")
        return False


def write_final_results(results, designs, output_dir, chai1_out):
    """Save final top-K results: rename structures and write CSV with sequences.

    For each top-K design:
      1. Copies the best-model CIF from Chai-1 output, converting to PDB
         with the naming scheme: rankNNN_{tag}_rf3_{orig_name}.pdb
      2. Reads chain sequences from the design's unique FASTA
      3. Writes everything to final_topk.csv with per-chain sequence columns

    Args:
        results: list of result dicts (sorted by combined score descending)
        designs: list of design dicts (from deduplicate_designs, with _tag set)
        output_dir: protocol output directory
        chai1_out: Chai-1 prediction output directory
    """
    from design_utils import read_fasta_chains

    top_k_dir = os.path.join(output_dir, 'topk_structures')
    os.makedirs(top_k_dir, exist_ok=True)

    # Lookup: design _tag → design dict
    design_by_tag = {d.get('_tag'): d for d in designs if d.get('_tag')}

    csv_path = os.path.join(output_dir, 'final_topk.csv')

    # Collect all chain IDs across top-K designs for CSV columns
    all_chain_ids = OrderedDict()
    for r in results:
        tag = r['design_name']
        d = design_by_tag.get(tag)
        if d and d.get('_unique_fasta') and os.path.exists(d['_unique_fasta']):
            chains = read_fasta_chains(d['_unique_fasta'])
            for cid in chains:
                all_chain_ids[cid] = None

    # Build fieldnames
    fieldnames = [
        'rank', 'design_name', 'orig_name',
        'best_model', 'plddt', 'ptm', 'iptm', 'combined',
        'structure_file',
    ]
    for cid in all_chain_ids:
        fieldnames.append(f'chain_{cid}')

    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for i, r in enumerate(results, 1):
            tag = r['design_name']
            orig_name = r.get('orig_name', '')
            best_model = r.get('best_model', '')

            # --- rename best-model structure ---
            out_name = f"rank{i:03d}_{tag}_rf3_{orig_name}.pdb"
            out_path = os.path.join(top_k_dir, out_name)

            if best_model != '' and not os.path.exists(out_path):
                cif_path = os.path.join(chai1_out, tag, f"pred.model_idx_{best_model}.cif")
                if os.path.exists(cif_path):
                    _cif_to_pdb(cif_path, out_path)

            # --- read chain sequences from FASTA ---
            chain_seqs = {}
            d = design_by_tag.get(tag)
            if d and d.get('_unique_fasta') and os.path.exists(d['_unique_fasta']):
                chain_seqs = read_fasta_chains(d['_unique_fasta'])

            # --- build CSV row ---
            row = {
                'rank': i,
                'design_name': tag,
                'orig_name': orig_name,
                'best_model': best_model,
                'plddt': f"{r['plddt']:.3f}" if r.get('plddt') else '',
                'ptm': f"{r['ptm']:.4f}" if r.get('ptm') else '',
                'iptm': f"{r['iptm']:.4f}" if r.get('iptm') else '',
                'combined': f"{r['combined']:.3f}" if r.get('combined') else '',
                'structure_file': out_name if os.path.exists(out_path) else '',
            }
            for cid in all_chain_ids:
                row[f'chain_{cid}'] = chain_seqs.get(cid, '')

            writer.writerow(row)

    n = len(results)
    print(f'\n[Phase 4] Top {n} structures saved to {top_k_dir}/')
    print(f'  Results saved to {csv_path}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    cfg = load_config(args.config)
    os.makedirs(args.output, exist_ok=True)

    # Phase 1: RFDiffusion3 backbone generation (N rounds, parallel)
    round_dirs = []
    rfd3_needed = []
    for r in range(1, args.num_rounds + 1):
        rfd3_out = os.path.join(args.output, f'rfd3_round_{r}')
        scores_csv = os.path.join(rfd3_out, 'rfd3_scores.csv')
        if args.skip_existing and os.path.exists(scores_csv):
            print(f'[Phase 1] Skipping round {r}/{args.num_rounds} — already completed')
            round_dirs.append(rfd3_out)
        else:
            rfd3_needed.append((r, rfd3_out))

    if rfd3_needed:
        processes = []
        for r, rfd3_out in rfd3_needed:
            print(f'[Phase 1] Submitting RFDiffusion3 round {r}/{args.num_rounds} '
                  f'({args.num_designs} designs)')
            cmd = build_rfd3_cmd(args, rfd3_out, args.config)
            proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
            processes.append((r, rfd3_out, proc))

        print(f'\n  Submitted {len(processes)} round(s) — waiting for all to complete...\n')
        for r, rfd3_out, proc in processes:
            stdout, stderr = proc.communicate()
            if proc.returncode == 0:
                if rfd3_out not in round_dirs:
                    round_dirs.append(rfd3_out)
                print(f'  Round {r} completed successfully.')
            else:
                print(f'  WARNING: round {r} failed (return code {proc.returncode})')
                if stderr:
                    for line in stderr.strip().splitlines()[-3:]:
                        print(f'    {line}')

    if not round_dirs:
        print('ERROR: all RFD3 rounds failed.')
        sys.exit(1)

    # Phase 2: Collect all designs, deduplicate, select top M
    all_designs, chain_meta = collect_designs_from_rounds(round_dirs)
    if not all_designs:
        print('ERROR: no designs collected from any round.')
        sys.exit(1)

    unique = deduplicate_designs(all_designs, chain_meta)

    # Select top M unique designs
    top_m = unique[:args.top_m]
    print(f'\n[Phase 2] Selected top {len(top_m)}/{len(unique)} unique designs:')
    for i, d in enumerate(top_m):
        info = f'  {i+1}. {d["design_name"]}  score={d["_score"]:.1f}'
        if d.get('_dedup_group_size', 1) > 1:
            info += f'  (from {d["_dedup_group_size"]} duplicates)'
        print(info)

    if not top_m:
        print('No designs selected. Exiting.')
        return

    # Write unique FASTA files for Chai-1
    unique_fasta_dir = write_unique_fastas(top_m,
                                           os.path.join(args.output, 'unique_designs'))

    # Phase 3: Chai-1 prediction
    results = run_chai1_batch(top_m, unique_fasta_dir,
                              os.path.join(args.output, 'chai1_preds'),
                              cfg, args)

    # Phase 4: Final ranking
    results.sort(key=lambda r: r['combined'], reverse=True)
    top_k = results[:args.top_k]

    print(f'\n[Phase 4] Top {len(top_k)} final designs:')
    for i, r in enumerate(top_k):
        print(f'  {i+1}. {r["design_name"]:30s} '
              f'pLDDT={r["plddt"]:.3f}  pTM={r["ptm"]:.3f}  iPTM={r["iptm"]:.3f}  '
              f'combined={r["combined"]:.3f}')

    chai1_out = os.path.join(args.output, 'chai1_preds')
    write_final_results(top_k, top_m, args.output, chai1_out)
    print('\nDone.')


if __name__ == '__main__':
    main()
