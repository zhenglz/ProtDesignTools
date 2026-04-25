#!/usr/bin/env python3
"""
Protocol: RFDiffusion3 → Chai-1 → ProteinMPNN → ESM-IF

Pipeline:
  1. RFDiffusion3: backbone generation for fragment design
  2. Chai-1: structure prediction (SLURM), select top M
  3. ProteinMPNN: for each top structure, L iterations of sequence design
     with 1-k randomly selected residues from the designed fragment
  4. ESM-IF: score designed sequences and rank
  5. Report best sequences in CSV

Usage:
  python protocols/fragment_design_mpnn.py \\
    --pdb input.pdb --design-regions "A90-120:20-30" \\
    --num-designs 50 --chai1-top 5 --mpnn-k 5 --mpnn-l 10 \\
    --output ./results
"""

import argparse, csv, os, random, re, sys, subprocess as sp, time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / 'tools'))

from _config import load_config
from chai1_tool import submit_chai1_slurm, check_job_status, extract_all_scores
from proteinmpnn_tool import run_design, parse_design_output, pdb2fasta

RFD3_SCRIPT = str(REPO_ROOT / 'tools' / 'rfdiffusion3_tool.py')

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def composite_rank(d):
    return (int(d.get('n_chainbreaks', 0)) * 100 +
            int(d.get('n_clashing', 0)) * 10 +
            abs(float(d.get('radius_of_gyration', 0)) - 12) * 0.1 +
            float(d.get('alanine_content', 0)) * 10 +
            abs(float(d.get('glycine_content', 0)) - 0.05) * 10 +
            float(d.get('max_ca_deviation', 0)))


def parse_design_regions(spec):
    """Parse 'A90-120:20-30' into [(chain, start, end), ...]."""
    regions = []
    for part in spec.split(';'):
        m = re.match(r'^([A-Za-z])(\d+)-(\d+):(\d+)-(\d+)$', part.strip())
        if not m:
            raise ValueError(f"Cannot parse region: '{part}'")
        regions.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return regions


def random_positions(regions, k):
    """Randomly select 1-k residues from the designed regions.

    Returns a positions string like "90A,92A-93A,95A" for ProteinMPNN.
    """
    all_res = []
    for chain, start, end in regions:
        all_res.extend([(chain, r) for r in range(start, end + 1)])
    if not all_res:
        return ''
    n = random.randint(1, min(k, len(all_res)))
    selected = sorted(random.sample(all_res, n), key=lambda x: (x[0], x[1]))
    # Build compact representation with ranges
    parts, cur_chain, cur_start, cur_end = [], None, None, None
    for chain, resi in selected:
        if cur_chain == chain and cur_end is not None and resi == cur_end + 1:
            cur_end = resi
        else:
            if cur_chain is not None:
                parts.append(f"{cur_start}{cur_chain}-{cur_end}{cur_chain}" if cur_start != cur_end
                             else f"{cur_start}{cur_chain}")
            cur_chain, cur_start, cur_end = chain, resi, resi
    if cur_chain is not None:
        parts.append(f"{cur_start}{cur_chain}-{cur_end}{cur_chain}" if cur_start != cur_end
                     else f"{cur_start}{cur_chain}")
    return ','.join(parts)


def find_mutations(wt_seq, mut_seq):
    """Compare two sequences, return list of mutation strings like 'M198F'."""
    if len(wt_seq) != len(mut_seq):
        return []
    return [f"{w}{i+1}{m}" for i, (w, m) in enumerate(zip(wt_seq, mut_seq)) if w != m]


def cif_to_pdb(cif_path, pdb_path):
    """Convert CIF to PDB using BioPython."""
    try:
        from Bio.PDB import MMCIFParser, PDBIO
        parser = MMCIFParser(QUIET=True)
        s = parser.get_structure('m', cif_path)
        io = PDBIO(); io.set_structure(s); io.save(pdb_path)
        return True
    except Exception as e:
        print(f"    CIF→PDB failed: {e}")
        return False


def esmif_score(pdb_file, mutations, esmif_python, esmif_script, work_dir):
    """Score list of mutations with ESM-IF. Returns list of (mutation, score)."""
    results = []
    for mut in mutations:
        cmd = [esmif_python, esmif_script, pdb_file, mut, work_dir]
        r = sp.run(cmd, capture_output=True, text=True)
        try:
            score = float(r.stdout.strip().split(': ')[-1])
            results.append((mut, score))
        except (ValueError, IndexError):
            print(f"    ESM-IF parse error for {mut}: {r.stdout[:100]}")
            results.append((mut, float('nan')))
    return results


# ---------------------------------------------------------------------------
# Phase helpers
# ---------------------------------------------------------------------------

def run_rfd3(args, rfd3_out, config_path):
    """Run RFDiffusion3 as subprocess."""
    cmd = [sys.executable, RFD3_SCRIPT, '--pdb', args.pdb,
           '--design-regions', args.design_regions,
           '--num-designs', str(args.num_designs),
           '--output', rfd3_out, '--top-n', '0',
           '--max-jobs', str(args.max_jobs)]
    if config_path:
        cmd += ['--config', config_path]
    if args.local:
        cmd.append('--local')
    else:
        cmd += ['--slurm-partition', args.slurm_partition, '--ncpus', str(args.ncpus)]
    print(f'[Phase 1] RFDiffusion3: {args.num_designs} designs')
    sp.run(cmd, check=True)


def run_chai1_batch(designs, rfd3_out, chai1_out, cfg, args):
    """Submit Chai-1 jobs for each design, wait, extract scores. Returns list of dicts."""
    os.makedirs(chai1_out, exist_ok=True)
    chai1_run = os.path.join(cfg['chai1']['chai1_dir'], 'run.sh')
    slurm_submit = cfg['slurm']['submit_script']
    fasta_dir = os.path.join(rfd3_out, 'top_designs_fastas')

    print(f'[Phase 2] Chai-1: {len(designs)} jobs')
    jobs = {}
    for d in designs:
        name = d['design_name']
        fasta = os.path.join(fasta_dir, f'{name}.fasta')
        if not os.path.exists(fasta):
            print(f'  SKIP {name}: no FASTA'); continue
        out_dir = os.path.join(chai1_out, name)
        jid = submit_chai1_slurm(fasta, out_dir, slurm_partition=args.slurm_partition,
                                 ncpus=args.ncpus, chai1_run=chai1_run,
                                 slurm_submit=slurm_submit)
        if jid:
            jobs[jid] = (name, out_dir)

    total = len(jobs)
    while jobs:
        done = []
        for jid, (name, out_dir) in jobs.items():
            st = check_job_status(jid)
            if st in ('COMPLETED', 'CD', 'COMPLETING'):
                done.append(jid)
            elif st in ('FAILED', 'F', 'CANCELLED', 'CA', 'TIMEOUT', 'TO'):
                print(f'  {name}: FAILED ({st})'); done.append(jid)
        for jid in done:
            del jobs[jid]
        if jobs:
            print(f'  [{time.strftime("%H:%M:%S")}] {total - len(jobs)}/{total} done, '
                  f'{len(jobs)} pending')
            time.sleep(30)

    results = []
    for d in designs:
        name, out_dir = d['design_name'], os.path.join(chai1_out, d['design_name'])
        scores = extract_all_scores(name, out_dir)
        if scores:
            results.append({'design_name': name, 'out_dir': out_dir,
                            'plddt': scores.get('best_plddt', 0),
                            'iptm': scores.get('best_iptm', 0),
                            'best_model': scores.get('best_model', 0)})
    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='RFD3 → Chai-1 → ProteinMPNN → ESM-IF',
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument('--pdb', required=True)
    p.add_argument('--design-regions', required=True,
                   help='Fragment spec: <chain><start>-<end>:<len_min>-<len_max>[;...]')
    p.add_argument('--num-designs', type=int, default=50,
                   help='RFD3 designs (default: 50)')
    p.add_argument('--chai1-top', type=int, default=5,
                   help='Top RFD3 designs to validate with Chai-1 (default: 5)')
    p.add_argument('--mpnn-k', type=int, default=5,
                   help='Max random residues per ProteinMPNN iteration (default: 5)')
    p.add_argument('--mpnn-l', type=int, default=10,
                   help='ProteinMPNN design iterations per structure (default: 10)')
    p.add_argument('--mpnn-n', type=int, default=100,
                   help='Sequences per ProteinMPNN run (default: 100)')
    p.add_argument('--output', '-o', default='./protocol_out')
    p.add_argument('--slurm-partition', default='4090')
    p.add_argument('--ncpus', type=int, default=4)
    p.add_argument('--max-jobs', type=int, default=4)
    p.add_argument('--local', action='store_true')
    p.add_argument('--config', help='Path to config JSON')
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    cfg = load_config(args.config)
    os.makedirs(args.output, exist_ok=True)
    regions = parse_design_regions(args.design_regions)

    # Phase 1: RFDiffusion3
    rfd3_out = os.path.join(args.output, 'rfd3_designs')
    run_rfd3(args, rfd3_out, args.config)

    # Phase 2: Select top designs by composite score
    csv_path = os.path.join(rfd3_out, 'rfd3_scores.csv')
    with open(csv_path) as f:
        designs = list(csv.DictReader(f))
    for d in designs:
        d['_score'] = composite_rank(d)
    designs.sort(key=lambda d: d['_score'])
    top_chai1 = designs[:args.chai1_top]
    print(f'[Phase 2] Selected top {len(top_chai1)}/{len(designs)} designs for Chai-1')

    # Phase 3: Chai-1 structure prediction
    chai1_out = os.path.join(args.output, 'chai1_preds')
    chai1_results = run_chai1_batch(top_chai1, rfd3_out, chai1_out, cfg, args)
    chai1_results.sort(key=lambda r: r['plddt'], reverse=True)

    if not chai1_results:
        print('No Chai-1 results. Exiting.'); return

    # Phase 4: ProteinMPNN + ESM-IF
    mpnn_out = os.path.join(args.output, 'mpnn_designs')
    os.makedirs(mpnn_out, exist_ok=True)
    esmif_work = os.path.join(args.output, 'esmif_scores')
    os.makedirs(esmif_work, exist_ok=True)

    # Load ProteinMPNN config
    pkg_dpath = cfg['proteinmpnn']['package_dpath']
    python_exe = f"{pkg_dpath}/../python_env/proteinmpnn/bin/python"

    # Load ESM-IF config
    esmif_python = cfg['esmif']['python_exe']
    esmif_script = cfg['esmif']['script_path']

    all_designs = []  # (chai1_name, iter, seq, mpnn_score, esmif_score, mutations)

    print(f'\n[Phase 4] ProteinMPNN + ESM-IF on top {len(chai1_results)} structures')
    for cr in chai1_results:
        name = cr['design_name']
        # Convert best CIF to PDB
        best_cif = os.path.join(cr['out_dir'], f"pred.model_idx_{cr['best_model']}.cif")
        pdb_file = os.path.join(mpnn_out, f"{name}.pdb")
        if not os.path.exists(pdb_file):
            if not cif_to_pdb(best_cif, pdb_file):
                print(f'  SKIP {name}: CIF→PDB conversion failed'); continue

        # Extract wild-type sequence from PDB
        wt_seq = pdb2fasta(pdb_file)
        if not wt_seq:
            print(f'  SKIP {name}: cannot extract PDB sequence'); continue

        print(f'\n  {name} (pLDDT={cr["plddt"]:.3f}, iPTM={cr["iptm"]:.3f})')

        for iteration in range(args.mpnn_l):
            pos_str = random_positions(regions, args.mpnn_k)
            if not pos_str:
                continue
            design_dir = os.path.join(mpnn_out, f"{name}_iter{iteration:02d}")
            os.makedirs(design_dir, exist_ok=True)

            success = run_design(pdb_file, pos_str, None, design_dir,
                                 num_seqs=args.mpnn_n, temperature=0.1,
                                 package_dpath=pkg_dpath, python_exe=python_exe)
            if not success:
                continue

            df = parse_design_output(design_dir, wt_seq)
            if df is None or len(df) < 2:
                continue

            # Take best designed sequence by ProteinMPNN score
            designs_df = df[df['name'] != 'wild_type'].sort_values('score')
            best = designs_df.iloc[0]
            design_seq = best['sequence']

            # Find mutations vs wild-type
            mutations = find_mutations(wt_seq, design_seq)
            if not mutations:
                continue

            # Score with ESM-IF
            esmif_results = esmif_score(pdb_file, mutations, esmif_python, esmif_script, esmif_work)
            esmif_total = sum(s for _, s in esmif_results if not (s != s))  # skip NaN
            n_scored = sum(1 for _, s in esmif_results if not (s != s))

            all_designs.append({
                'chai1_design': name,
                'iteration': iteration,
                'positions': pos_str,
                'n_mutations': len(mutations),
                'n_scored': n_scored,
                'mpnn_score': best['score'],
                'esmif_total': esmif_total,
                'esmif_mean': esmif_total / n_scored if n_scored else float('nan'),
                'sequence': design_seq,
                'mutations': ','.join(mutations),
                'plddt': cr['plddt'],
                'iptm': cr['iptm'],
            })
            print(f'    iter {iteration:02d}: positions={pos_str}, '
                  f'mutations={len(mutations)}, esmif_total={esmif_total:.3f}')

    # Phase 5: Report
    if not all_designs:
        print('\nNo designs generated. Exiting.'); return

    all_designs.sort(key=lambda d: d['esmif_total'])
    print(f'\n[Phase 5] Top designs by ESM-IF score (lower = better):')
    for i, d in enumerate(all_designs[:20]):
        print(f'  {i+1:3d}. {d["chai1_design"]:30s} iter={d["iteration"]:02d}  '
              f'mpnn={d["mpnn_score"]:.3f}  esmif_total={d["esmif_total"]:.3f}  '
              f'mutations={d["n_mutations"]}  pos={d["positions"]}')

    result_csv = os.path.join(args.output, 'final_designs.csv')
    with open(result_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=[
            'rank', 'chai1_design', 'iteration', 'positions', 'n_mutations',
            'n_scored', 'mpnn_score', 'esmif_total', 'esmif_mean',
            'sequence', 'mutations', 'plddt', 'iptm'])
        w.writeheader()
        for i, d in enumerate(all_designs, 1):
            d['rank'] = i
            w.writerow(d)
    print(f'\nResults saved to {result_csv}')


if __name__ == '__main__':
    main()
