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
from collections import OrderedDict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / 'tools'))

from _config import load_config
from chai1_tool import submit_chai1_slurm, check_job_status, extract_all_scores
from design_utils import (
    composite_rank, collect_designs_from_rounds,
    deduplicate_designs, write_unique_fastas,
    get_designed_chain_ids,
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
    p.add_argument('--weights', default='plddt:1;iptm:1;helicity:1',
                   help='Weighted ranking metrics and weights, semicolon-separated '
                        '(default: plddt:1;iptm:1;helicity:1)')
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
    n_total = len(designs)
    print(f'  Extracting scores for {n_total} design(s)...')
    results = []
    for idx, d in enumerate(designs, 1):
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
        if idx % 50 == 0 or idx == n_total:
            print(f'    [{time.strftime("%H:%M:%S")}] {idx}/{n_total} scores extracted')
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


def calc_helicity(cif_path, designed_chains, python_exe=None):
    """Calculate alpha-helix fraction for designed chains using mdtraj DSSP.

    Tries direct ``import mdtraj`` first.  If unavailable, falls back
    to calling *python_exe* (e.g. the sfct conda environment) as a
    subprocess.  As a last resort uses a BioPython CA-distance heuristic.

    Args:
        cif_path: Path to the CIF structure file.
        designed_chains: Set of chain IDs (e.g. {'C'}) to measure.
        python_exe: Optional path to a Python with mdtraj installed.

    Returns:
        Float in [0, 1], or 0.0 on total failure.
    """
    # --- attempt 1: direct import ---
    try:
        import mdtraj as md
        return _helicity_mdtraj(cif_path, designed_chains, md)
    except ImportError:
        pass

    # --- attempt 2: subprocess via python_exe (e.g. sfct env) ---
    if python_exe and os.path.exists(python_exe):
        try:
            return _helicity_subprocess(cif_path, designed_chains, python_exe)
        except Exception:
            pass

    # --- attempt 3: BioPython fallback ---
    try:
        from Bio.PDB import MMCIFParser
        return _helicity_biopython(cif_path, designed_chains, MMCIFParser)
    except ImportError:
        pass

    return 0.0


def _helicity_mdtraj(cif_path, designed_chains, md):
    """mdtraj DSSP-based helicity for designed chains."""
    traj = md.load(cif_path)
    dssp_per_res = md.compute_dssp(traj, simplified=True)[0]
    topology = traj.topology
    n_helix = n_total = 0
    for res in topology.residues:
        cid = str(res.chain.chain_id)
        if cid not in designed_chains:
            continue
        n_total += 1
        if dssp_per_res[res.index] == 'H':
            n_helix += 1
    return n_helix / n_total if n_total else 0.0


def _helicity_subprocess(cif_path, designed_chains, python_exe):
    """Run mdtraj DSSP in a subprocess using a dedicated python."""
    import json, textwrap
    code = textwrap.dedent(f"""\
        import mdtraj as md, json, sys
        traj = md.load({json.dumps(cif_path)})
        dssp = md.compute_dssp(traj, simplified=True)[0]
        top = traj.topology
        designed = {json.dumps(sorted(designed_chains))}
        n_helix = n_total = 0
        for res in top.residues:
            cid = str(res.chain.chain_id)
            if cid not in designed:
                continue
            n_total += 1
            if dssp[res.index] == 'H':
                n_helix += 1
        print(n_helix / n_total if n_total else 0.0)
    """)
    result = sp.run([python_exe, '-c', code],
                    stdout=sp.PIPE, stderr=sp.PIPE,
                    universal_newlines=True, timeout=120)
    if result.returncode == 0 and result.stdout.strip():
        return float(result.stdout.strip())
    raise RuntimeError(f"subprocess failed: {result.stderr}")


def _helicity_biopython(cif_path, designed_chains, MMCIFParser):
    """CA-CA[i+4] distance heuristic for designed chains."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("model", cif_path)
    n_helix = n_total = 0
    for model in structure:
        for chain in model:
            if chain.id not in designed_chains:
                continue
            residues = [r for r in chain if r.id[0] == ' ']
            for i, res in enumerate(residues):
                n_total += 1
                if i + 4 >= len(residues):
                    continue
                try:
                    ca_i = res['CA'].get_vector()
                    ca_i4 = residues[i + 4]['CA'].get_vector()
                    if (ca_i - ca_i4).norm() < 7.0:
                        n_helix += 1
                except (KeyError, IndexError):
                    continue
    return n_helix / n_total if n_total else 0.0


def write_final_results(all_results, designs, output_dir, chai1_out,
                        top_k, weights_str, designed_chains=None, python_exe=None):
    """Rank all results by weighted score, select top K, rename structures,
    write CSV with scores, helicity, and chain sequences.

    Args:
        all_results: list of result dicts from run_chai1_batch
        designs: list of design dicts (from deduplicate_designs, with _tag)
        output_dir: protocol output directory
        chai1_out: Chai-1 prediction output directory
        top_k: number of top designs to keep
        weights_str: weighting spec, e.g. "plddt:1;iptm:1;helicity:1"
    """
    from design_utils import read_fasta_chains

    # --- parse weights ---
    weights = {}
    for token in weights_str.replace(',', ';').split(';'):
        token = token.strip()
        if not token:
            continue
        parts = token.split(':')
        if len(parts) == 2:
            key, val = parts[0].strip().lower(), float(parts[1])
            if val != 0:
                weights[key] = val
    if not weights:
        weights = {'plddt': 1.0, 'iptm': 1.0, 'helicity': 1.0}
    norm = sum(weights.values())

    # --- calculate helicity and weighted score for every result ---
    for r in all_results:
        tag = r['design_name']
        # Read helicity from best-model CIF
        helicity = 0.0
        bm = r.get('best_model', '')
        if bm != '':
            cif_path = os.path.join(chai1_out, tag, f"pred.model_idx_{bm}.cif")
            if os.path.exists(cif_path):
                helicity = calc_helicity(cif_path, designed_chains or set(), python_exe)
        r['helicity'] = helicity

        # Weighted combined score (normalised)
        wscore = 0.0
        for key, w in weights.items():
            if key == 'plddt':
                wscore += w * r.get('plddt', 0)
            elif key == 'iptm':
                wscore += w * (r.get('iptm', 0) * 100)  # iptm 0-1 scale → percentage
            elif key == 'helicity':
                wscore += w * (helicity * 100)  # helicity 0-1 → percentage
            elif key == 'ptm':
                wscore += w * (r.get('ptm', 0) * 100)
        r['weighted_score'] = wscore / norm if norm else 0

    # --- rank by weighted score, select top K ---
    all_results.sort(key=lambda r: r['weighted_score'], reverse=True)
    top = all_results[:top_k]

    # --- build lookup ---
    design_by_tag = {d.get('_tag'): d for d in designs if d.get('_tag')}

    top_k_dir = os.path.join(output_dir, 'final_results', 'topk_structures')
    os.makedirs(top_k_dir, exist_ok=True)

    # Collect chain IDs across top designs for CSV columns
    all_chain_ids = OrderedDict()
    for r in top:
        tag = r['design_name']
        d = design_by_tag.get(tag)
        if d and d.get('_unique_fasta') and os.path.exists(d['_unique_fasta']):
            chains = read_fasta_chains(d['_unique_fasta'])
            for cid in chains:
                all_chain_ids[cid] = None

    fieldnames = [
        'rank', 'design_name', 'orig_name',
        'best_model', 'plddt', 'ptm', 'iptm', 'helicity',
        'weighted_score', 'structure_file',
    ]
    for cid in all_chain_ids:
        fieldnames.append(f'chain_{cid}')

    csv_path = os.path.join(output_dir, 'final_results', 'final_topk.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for i, r in enumerate(top, 1):
            tag = r['design_name']
            orig_name = r.get('orig_name', '')
            bm = r.get('best_model', '')

            # --- rename best-model structure ---
            out_name = f"rank{i:03d}_{tag}_rf3_{orig_name}.pdb"
            out_path = os.path.join(top_k_dir, out_name)
            if bm != '' and not os.path.exists(out_path):
                cif_path = os.path.join(chai1_out, tag, f"pred.model_idx_{bm}.cif")
                if os.path.exists(cif_path):
                    _cif_to_pdb(cif_path, out_path)

            # --- read chain sequences ---
            chain_seqs = {}
            d = design_by_tag.get(tag)
            if d and d.get('_unique_fasta') and os.path.exists(d['_unique_fasta']):
                chain_seqs = read_fasta_chains(d['_unique_fasta'])

            row = {
                'rank': i,
                'design_name': tag,
                'orig_name': orig_name,
                'best_model': bm,
                'plddt': f"{r['plddt']:.3f}" if r.get('plddt') else '',
                'ptm': f"{r['ptm']:.4f}" if r.get('ptm') else '',
                'iptm': f"{r['iptm']:.4f}" if r.get('iptm') else '',
                'helicity': f"{r['helicity']:.4f}",
                'weighted_score': f"{r['weighted_score']:.3f}",
                'structure_file': out_name if os.path.exists(out_path) else '',
            }
            for cid in all_chain_ids:
                row[f'chain_{cid}'] = chain_seqs.get(cid, '')
            writer.writerow(row)

    # --- print summary ---
    print(f'\n[Phase 4] Top {len(top)} final designs (weighted by {weights_str}):')
    for i, r in enumerate(top, 1):
        print(f'  {i}. {r["design_name"]:30s}  '
              f'pLDDT={r["plddt"]:.3f}  iPTM={r["iptm"]:.4f}  '
              f'helicity={r["helicity"]:.3f}  '
              f'score={r["weighted_score"]:.3f}')
    print(f'\n  Structures: {top_k_dir}/')
    print(f'  CSV: {csv_path}')


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
        rfd3_out = os.path.join(args.output, 'step1_rfd3', f'round_{r}')
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
                                           os.path.join(args.output, 'step1_rfd3', 'unique_designs'))

    # Phase 3: Chai-1 prediction
    results = run_chai1_batch(top_m, unique_fasta_dir,
                              os.path.join(args.output, 'step2_chai1'),
                              cfg, args)

    # Phase 4: Final ranking (weighted score with helicity)
    chai1_out = os.path.join(args.output, 'step2_chai1')
    designed_chains = get_designed_chain_ids(chain_meta)
    python_exe = cfg.get('esmif', {}).get('python_executable')  # sfct env with mdtraj
    write_final_results(results, top_m, args.output, chai1_out,
                        args.top_k, args.weights,
                        designed_chains=designed_chains, python_exe=python_exe)
    print('\nDone.')


if __name__ == '__main__':
    main()
