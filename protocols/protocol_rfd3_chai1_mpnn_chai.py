#!/usr/bin/env python3
"""
Protocol: RFDiffusion3 → Chai-1 → ProteinMPNN → Chai-1 validation

Pipeline:
  1. RFDiffusion3: backbone generation (N rounds)
  2. Chai-1: predict structures for top unique designs, rank
  3. ProteinMPNN: redesign designed/unfixed regions on top-K structures
  4. Chai-1: predict structures for top MPNN designs, final ranking

Usage:
  python protocols/protocol_rfd3_chai1_mpnn_chai.py \\
    --pdb input.pdb --design-regions "C820-831:5-15" \\
    --num-rounds 3 --num-designs 50 --top-m 20 --top-k 3 \\
    --mpnn-num-seqs 100 --mpnn-top-n 10 --output ./results
"""

import argparse, csv, os, re, sys, subprocess as sp, time, shutil
from collections import OrderedDict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / 'tools'))

from _config import load_config
from chai1_tool import submit_chai1_slurm, check_job_status, extract_all_scores
from design_utils import (
    composite_rank, collect_designs_from_rounds,
    deduplicate_designs, write_unique_fastas,
    get_designed_chain_ids, read_fasta_chains,
)
from proteinmpnn_tool import run_design, parse_design_output

RFD3_SCRIPT = str(REPO_ROOT / 'tools' / 'rfdiffusion3_tool.py')

# Standard 3-letter to 1-letter amino acid codes
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


# ---------------------------------------------------------------------------
# Shared helpers (from design_validate.py)
# ---------------------------------------------------------------------------

def _cif_to_pdb(cif_path, pdb_path):
    """Convert CIF to PDB using BioPython."""
    try:
        from Bio.PDB import MMCIFParser, PDBIO
        parser = MMCIFParser(QUIET=True)
        s = parser.get_structure('m', cif_path)
        io = PDBIO()
        io.set_structure(s)
        io.save(pdb_path)
        return True
    except ImportError:
        shutil.copy2(cif_path, pdb_path)
        return False
    except Exception as e:
        print(f"    CIF->PDB failed: {e}")
        return False


def calc_helicity(cif_path, designed_chains, python_exe=None):
    """Calculate alpha-helix fraction for designed chains using mdtraj DSSP."""
    # --- attempt 1: direct import ---
    try:
        import mdtraj as md
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
    except ImportError:
        pass

    # --- attempt 2: subprocess via python_exe ---
    if python_exe and os.path.exists(python_exe):
        try:
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
        except Exception:
            pass

    # --- attempt 3: BioPython fallback ---
    try:
        from Bio.PDB import MMCIFParser
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
    except ImportError:
        pass

    return 0.0


def parse_weights(weights_str):
    """Parse 'plddt:1;iptm:1;helicity:1' into {key: weight, ...} dict."""
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
    return weights


# ---------------------------------------------------------------------------
# Design-region helpers
# ---------------------------------------------------------------------------

def parse_design_regions(spec):
    """Parse 'C820-831:5-15' into [(chain, start, end), ...]."""
    regions = []
    for part in spec.split(';'):
        m = re.match(r'^([A-Za-z])(\d+)-(\d+):(\d+)-(\d+)$', part.strip())
        if not m:
            raise ValueError(f"Cannot parse region: '{part}'")
        regions.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return regions


def get_chain_seq(pdb_file, chain_id):
    """Extract 1-letter amino acid sequence for a chain from PDB."""
    try:
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('x', pdb_file)
        for chain in structure[0]:
            if chain.id == chain_id:
                residues = [r for r in chain
                            if r.id[0] == ' ' and r.resname in THREE_TO_ONE]
                return ''.join(THREE_TO_ONE[r.resname] for r in residues)
    except Exception:
        pass
    return ''


def get_all_chain_seqs(pdb_file):
    """Extract sequences for all chains in a PDB file. Returns OrderedDict."""
    seqs = OrderedDict()
    try:
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('x', pdb_file)
        for chain in structure[0]:
            residues = [r for r in chain
                        if r.id[0] == ' ' and r.resname in THREE_TO_ONE]
            if residues:
                seqs[chain.id] = ''.join(THREE_TO_ONE[r.resname] for r in residues)
    except Exception:
        pass
    return seqs


def design_regions_to_mpnn_positions(regions):
    """Convert design regions to MPNN position string.
    [(chain, start, end), ...] -> "820C-831C"
    """
    parts = []
    for chain, start, end in regions:
        parts.append(f"{start}{chain}-{end}{chain}")
    return ','.join(parts)


def find_mutations(wt_seq, mut_seq):
    """Compare two sequences, return list of mutation strings like 'M198F'."""
    if len(wt_seq) != len(mut_seq):
        return []
    return [f"{w}{i+1}{m}" for i, (w, m) in enumerate(zip(wt_seq, mut_seq)) if w != m]


def apply_mutations_to_chains(chain_seqs, mutations):
    """Apply mutation strings to per-chain sequences.

    Args:
        chain_seqs: OrderedDict of chain_id -> sequence (sorted order = concatenation order)
        mutations: list of mutation strings like 'M198F' (1-indexed in concatenated sequence)

    Returns:
        OrderedDict with mutated sequences
    """
    # Build offset map
    offsets = {}
    offset = 0
    for ch, seq in chain_seqs.items():
        offsets[ch] = offset
        offset += len(seq)
    total_len = offset

    new_seqs = OrderedDict(chain_seqs)
    for m_str in mutations:
        if len(m_str) < 3:
            continue
        pos_1based = int(m_str[1:-1])  # 1-indexed in concatenated seq
        new_aa = m_str[-1]
        if pos_1based < 1 or pos_1based > total_len:
            continue
        pos = pos_1based - 1  # 0-indexed

        # Find which chain this position falls in
        for ch in new_seqs:
            ch_len = len(new_seqs[ch])
            if offsets[ch] <= pos < offsets[ch] + ch_len:
                rel = pos - offsets[ch]
                seq_list = list(new_seqs[ch])
                seq_list[rel] = new_aa
                new_seqs[ch] = ''.join(seq_list)
                break
    return new_seqs


def build_chai1_fasta(chain_seqs):
    """Build Chai-1 format FASTA string from chain_id -> sequence dict."""
    lines = []
    for cid, seq in chain_seqs.items():
        lines.append(f">protein|{cid}")
        lines.append(seq)
    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# Phase 1: RFDiffusion3
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


# ---------------------------------------------------------------------------
# Phase 2: Chai-1 batch
# ---------------------------------------------------------------------------

def run_chai1_batch(designs, fasta_dir, chai1_out, cfg, args):
    """Submit Chai-1 SLURM jobs for each design, wait, extract scores.

    Returns list of result dicts.
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
    jobs = {}
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

    # Wait for jobs
    total = len(jobs)
    while jobs:
        done = []
        for jid, (tag, out_dir) in jobs.items():
            still_running = check_job_status(jid)
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


# ---------------------------------------------------------------------------
# Phase 2 ranking (step 2 results -> final_results/step2)
# ---------------------------------------------------------------------------

def rank_and_save_step2(results, designs, output_dir, chai1_out,
                        top_k, weights_str, designed_chains=None, python_exe=None):
    """Rank step 2 Chai-1 results, write to final_results/step2/."""
    weights = parse_weights(weights_str)
    norm = sum(weights.values())

    # Calculate helicity and weighted score
    for r in results:
        tag = r['design_name']
        helicity = 0.0
        bm = r.get('best_model', '')
        if bm != '':
            cif_path = os.path.join(chai1_out, tag, f"pred.model_idx_{bm}.cif")
            if os.path.exists(cif_path):
                helicity = calc_helicity(cif_path, designed_chains or set(), python_exe)
        r['helicity'] = helicity

        wscore = 0.0
        for key, w in weights.items():
            if key == 'plddt':
                wscore += w * r.get('plddt', 0)
            elif key == 'iptm':
                wscore += w * (r.get('iptm', 0) * 100)
            elif key == 'helicity':
                wscore += w * (helicity * 100)
            elif key == 'ptm':
                wscore += w * (r.get('ptm', 0) * 100)
        r['weighted_score'] = wscore / norm if norm else 0

    # Rank by weighted score
    results.sort(key=lambda r: r['weighted_score'], reverse=True)
    top = results[:top_k]

    design_by_tag = {d.get('_tag'): d for d in designs if d.get('_tag')}

    # Write to final_results/step2/
    out_subdir = os.path.join(output_dir, 'final_results', 'step2')
    top_k_dir = os.path.join(out_subdir, 'topk_structures')
    os.makedirs(top_k_dir, exist_ok=True)

    # Collect chain IDs
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

    csv_path = os.path.join(out_subdir, 'final_topk.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for i, r in enumerate(top, 1):
            tag = r['design_name']
            orig_name = r.get('orig_name', '')
            bm = r.get('best_model', '')

            out_name = f"rank{i:03d}_{tag}_rf3_{orig_name}.pdb"
            out_path = os.path.join(top_k_dir, out_name)
            if bm != '' and not os.path.exists(out_path):
                cif_path = os.path.join(chai1_out, tag, f"pred.model_idx_{bm}.cif")
                if os.path.exists(cif_path):
                    _cif_to_pdb(cif_path, out_path)

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

    print(f'\n[Step 2] Top {len(top)} designs (weighted by {weights_str}):')
    for i, r in enumerate(top, 1):
        print(f'  {i}. {r["design_name"]:30s}  '
              f'pLDDT={r["plddt"]:.3f}  iPTM={r["iptm"]:.4f}  '
              f'helicity={r["helicity"]:.3f}  '
              f'score={r["weighted_score"]:.3f}')
    print(f'  Structures: {top_k_dir}/')
    print(f'  CSV: {csv_path}')

    return top  # return ranked top designs for Step 3


# ---------------------------------------------------------------------------
# Phase 3: ProteinMPNN design
# ---------------------------------------------------------------------------

def run_mpnn_phase(top_structures, output_dir, cfg, args, design_regions):
    """Run ProteinMPNN on top structures from step 2.

    Args:
        top_structures: list of ranked result dicts from rank_and_save_step2
        output_dir: protocol output directory
        cfg: config dict
        args: CLI args
        design_regions: list of (chain, start, end) tuples

    Returns:
        List of dicts with MPNN design results (tag, mpnn_tag, fasta_content, mpnn_score, ...)
    """
    mpnn_out = os.path.join(output_dir, 'step3_mpnn')
    os.makedirs(mpnn_out, exist_ok=True)

    # MPNN config
    pkg_dpath = cfg['proteinmpnn']['package_dpath']
    mpnn_python = f"{pkg_dpath}/../python_env/proteinmpnn/bin/python"

    # Designed chains
    designed_chains_set = sorted(set(r[0] for r in design_regions))

    # MPNN position string
    mpnn_positions = design_regions_to_mpnn_positions(design_regions)
    print(f'  MPNN design positions: {mpnn_positions}')

    all_mpnn_designs = []

    for rank_idx, struct in enumerate(top_structures, 1):
        tag = struct['design_name']
        bm = struct.get('best_model', '')
        chai1_out_in = os.path.join(output_dir, 'step2_chai1')
        cif_path = os.path.join(chai1_out_in, tag, f"pred.model_idx_{bm}.cif")

        if not os.path.exists(cif_path):
            print(f'  SKIP {tag}: best-model CIF not found')
            continue

        # Convert best CIF to PDB for MPNN
        pdb_file = os.path.join(mpnn_out, f"{tag}.pdb")
        if not os.path.exists(pdb_file):
            if not _cif_to_pdb(cif_path, pdb_file):
                print(f'  SKIP {tag}: CIF->PDB failed')
                continue

        # Get original sequences for all chains (needed for Chai-1 FASTA)
        all_chain_seqs = get_all_chain_seqs(pdb_file)
        if not all_chain_seqs:
            print(f'  SKIP {tag}: no chains extracted from PDB')
            continue

        # Get wild-type concatenated sequence for designed chains
        wt_parts = []
        for ch in designed_chains_set:
            if ch in all_chain_seqs:
                wt_parts.append(all_chain_seqs[ch])
        if not wt_parts:
            print(f'  SKIP {tag}: no designed-chain sequences found')
            continue
        wt_seq = ''.join(wt_parts)

        # Build MPNN exclude list: if not designing unfixed, keep all non-region positions fixed
        exclude_str = None
        if not args.mpnn_design_unfixed:
            # Calculate all residue positions NOT in design regions
            all_residues = []
            with open(pdb_file) as fh:
                for line in fh:
                    if line.startswith('ATOM'):
                        res_idx = int(line[22:26].strip())
                        chain_id = line[21]
                        all_residues.append((res_idx, chain_id))
            # Deduplicate
            all_residues = list(OrderedDict.fromkeys(all_residues))

            # Build set of designed positions
            designed_positions = set()
            for chain, start, end in design_regions:
                for r in range(start, end + 1):
                    designed_positions.add((r, chain))

            # Fixed = all residues not in designed positions
            fixed_parts = []
            for res_idx, chain_id in all_residues:
                if (res_idx, chain_id) not in designed_positions:
                    fixed_parts.append(f"{res_idx}{chain_id}")
            exclude_str = ','.join(fixed_parts) if fixed_parts else None
        else:
            # Designing all residues in the PDB — no exclude needed
            pass

        print(f'\n  [{rank_idx}] {tag}: running ProteinMPNN (pLDDT={struct.get("plddt",0):.3f})')
        if exclude_str:
            print(f'    Fixed positions: {len(exclude_str.split(","))} residues')

        # Run MPNN design
        design_dir = os.path.join(mpnn_out, tag)
        success = run_design(
            pdb_file, mpnn_positions, exclude_str, design_dir,
            num_seqs=args.mpnn_num_seqs, temperature=args.mpnn_temperature,
            package_dpath=pkg_dpath, python_exe=mpnn_python,
        )
        if not success:
            print(f'    SKIP {tag}: MPNN design failed')
            continue

        # Parse MPNN results
        df = parse_design_output(design_dir, wt_seq)
        if df is None or len(df) < 2:
            print(f'    SKIP {tag}: no MPNN results parsed')
            continue

        # Rank designed sequences by MPNN score (lower is better)
        designs_df = df[df['name'] != 'wild_type'].sort_values('score')
        n_designs = len(designs_df)
        print(f'    MPNN: {n_designs} designs, top score={designs_df.iloc[0]["score"]:.3f}')

        # Store all designs
        # Build full concatenated wild-type from ALL chains (MPNN outputs all chains)
        full_wt_seq = ''.join(all_chain_seqs.values())
        chain_offsets = {}
        off = 0
        for ch, seq in all_chain_seqs.items():
            chain_offsets[ch] = off
            off += len(seq)
        full_wt_len = off

        for mpnn_idx, (_, row) in enumerate(designs_df.iterrows()):
            mpnn_tag = f"{tag}_mpnn{mpnn_idx:03d}"
            design_seq_raw = row['sequence']

            # Clean separators that MPNN may insert between chains
            design_seq_full = design_seq_raw.replace(':', '').replace('/', '')
            design_len = len(design_seq_full)

            if design_len == full_wt_len:
                # Full multi-chain output — find all mutations by position
                new_seqs = OrderedDict(all_chain_seqs)
                n_mut = 0
                for ch in new_seqs:
                    ch_len = len(new_seqs[ch])
                    ch_off = chain_offsets[ch]
                    seq_list = list(new_seqs[ch])
                    for rel_pos in range(ch_len):
                        fp = ch_off + rel_pos
                        if fp < design_len and design_seq_full[fp] != full_wt_seq[fp]:
                            seq_list[rel_pos] = design_seq_full[fp]
                            n_mut += 1
                    new_seqs[ch] = ''.join(seq_list)
                mutations_found = n_mut > 0
                full_fasta_seqs = new_seqs
            else:
                # Only designed-chain output — use mutation approach
                mutations = find_mutations(wt_seq, design_seq_raw)
                if not mutations and mpnn_idx > 0:
                    continue  # skip identical sequences
                mutations_found = bool(mutations)
                new_chain_seqs = apply_mutations_to_chains(
                    OrderedDict((ch, all_chain_seqs[ch]) for ch in designed_chains_set
                                if ch in all_chain_seqs),
                    mutations,
                )
                full_fasta_seqs = OrderedDict()
                for ch, seq in all_chain_seqs.items():
                    full_fasta_seqs[ch] = new_chain_seqs.get(ch, seq)
                mutations_found = bool(mutations)

            if not mutations_found and mpnn_idx > 0:
                continue  # skip identical sequences (keep first = best)

            fasta_content = build_chai1_fasta(full_fasta_seqs)

            all_mpnn_designs.append({
                'mpnn_tag': mpnn_tag,
                'source_tag': tag,
                'rank': rank_idx,
                'mpnn_score': row['score'],
                'mpnn_recovery': row['recovery'],
                'design_seq': design_seq_raw,
                'n_mutations': sum(1 for a, b in zip(full_wt_seq, design_seq_full) if a != b)
                               if design_len == full_wt_len else 0,
                'fasta_content': fasta_content,
                'plddt': struct.get('plddt', 0),
                'iptm': struct.get('iptm', 0),
                'full_fasta_seqs': full_fasta_seqs,
            })

    if not all_mpnn_designs:
        print('\n[Step 3] No MPNN designs generated.')
        return []

    # Sort by MPNN score (lower = better)
    all_mpnn_designs.sort(key=lambda d: d['mpnn_score'])

    # Save all MPNN results
    csv_path = os.path.join(mpnn_out, 'mpnn_all_results.csv')
    with open(csv_path, 'w', newline='') as f:
        fieldnames = [
            'rank', 'mpnn_tag', 'source_tag', 'structure_rank',
            'mpnn_score', 'mpnn_recovery', 'n_mutations', 'design_seq',
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i, d in enumerate(all_mpnn_designs, 1):
            writer.writerow({
                'rank': i,
                'mpnn_tag': d['mpnn_tag'],
                'source_tag': d['source_tag'],
                'structure_rank': d['rank'],
                'mpnn_score': f"{d['mpnn_score']:.4f}",
                'mpnn_recovery': f"{d['mpnn_recovery']:.4f}",
                'n_mutations': d['n_mutations'],
                'design_seq': d['design_seq'],
            })

    print(f'\n[Step 3] Total MPNN designs: {len(all_mpnn_designs)}')
    print(f'  Results: {csv_path}')

    # Select top half (or --mpnn-chai1-top) for step 4 Chai-1 validation
    n_for_chai1 = args.mpnn_chai1_top if args.mpnn_chai1_top > 0 else max(1, len(all_mpnn_designs) // 2)
    top_mpnn = all_mpnn_designs[:n_for_chai1]

    # Write Chai-1 FASTA files for top MPNN designs
    fasta_dir = os.path.join(mpnn_out, 'chai1_fastas')
    os.makedirs(fasta_dir, exist_ok=True)
    for d in top_mpnn:
        fasta_path = os.path.join(fasta_dir, f"{d['mpnn_tag']}.fasta")
        with open(fasta_path, 'w') as f:
            f.write(d['fasta_content'] + '\n')
        d['_fasta_path'] = fasta_path

    print(f'  Submitted to Step 4 Chai-1: {len(top_mpnn)} designs')
    print(f'  FASTA files: {fasta_dir}/')

    return top_mpnn


# ---------------------------------------------------------------------------
# Phase 4: Chai-1 validation of MPNN designs
# ---------------------------------------------------------------------------

def run_mpnn_chai1(mpnn_designs, output_dir, cfg, args, designed_chains):
    """Run Chai-1 on MPNN-designed sequences.

    Args:
        mpnn_designs: list of MPNN design dicts from run_mpnn_phase
        output_dir: protocol output directory
        cfg: config dict
        args: CLI args
        designed_chains: set of chain IDs for helicity calculation

    Returns:
        List of ranked result dicts for final output.
    """
    if not mpnn_designs:
        return []

    chai1_out = os.path.join(output_dir, 'step4_chai1')
    os.makedirs(chai1_out, exist_ok=True)
    chai1_run = os.path.join(cfg['chai1']['chai1_dir'], 'run.sh')
    slurm_submit = cfg['slurm']['submit_script']
    python_exe = cfg.get('esmif', {}).get('python_executable')

    # Submit jobs
    to_submit = []
    skipped = []
    for d in mpnn_designs:
        tag = d['mpnn_tag']
        out_dir = os.path.join(chai1_out, tag)
        all_five = all(
            os.path.exists(os.path.join(out_dir, f'pred.model_idx_{m}.cif'))
            for m in range(5)
        )
        if args.skip_existing and all_five:
            skipped.append((tag, out_dir, d))
        else:
            to_submit.append(d)

    if skipped:
        print(f'  Skipping {len(skipped)} MPNN design(s) with existing Chai-1 results')

    jobs = {}
    if to_submit:
        print(f'[{time.strftime("%H:%M:%S")}] Submitting {len(to_submit)} Chai-1 job(s)...')
        for d in to_submit:
            tag = d['mpnn_tag']
            fasta = d.get('_fasta_path', '')
            if not fasta or not os.path.exists(fasta):
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
        print('No Chai-1 jobs submitted and no existing results.')
        return []

    # Wait for jobs
    total = len(jobs)
    while jobs:
        done = []
        for jid, (tag, out_dir) in jobs.items():
            still_running = check_job_status(jid)
            if not still_running:
                done.append(jid)
        for jid in done:
            del jobs[jid]
        if jobs:
            print(f'  [{time.strftime("%H:%M:%S")}] {total - len(jobs)}/{total} done, '
                  f'{len(jobs)} pending')
            time.sleep(30)

    # Extract scores for all designs (submitted + skipped)
    all_results = []
    for d in mpnn_designs:
        tag = d['mpnn_tag']
        out_dir = os.path.join(chai1_out, tag)
        scores = extract_all_scores(tag, out_dir)
        if scores:
            all_results.append({
                'design_name': tag,
                'source_tag': d['source_tag'],
                'mpnn_score': d['mpnn_score'],
                'n_mutations': d['n_mutations'],
                'design_seq': d['design_seq'],
                'best_model': scores.get('best_model', ''),
                'plddt': scores.get('best_plddt', 0),
                'ptm': scores.get('best_ptm', 0),
                'iptm': scores.get('best_iptm', 0),
                'combined': scores.get('best_plddt', 0) * scores.get('best_iptm', 0),
            })

    if not all_results:
        print('[Step 4] No Chai-1 results extracted.')
        return []

    # Calculate helicity and weighted score
    weights = parse_weights(args.weights)
    norm = sum(weights.values())

    for r in all_results:
        tag = r['design_name']
        helicity = 0.0
        bm = r.get('best_model', '')
        if bm != '':
            cif_path = os.path.join(chai1_out, tag, f"pred.model_idx_{bm}.cif")
            if os.path.exists(cif_path):
                helicity = calc_helicity(cif_path, designed_chains or set(), python_exe)
        r['helicity'] = helicity

        wscore = 0.0
        for key, w in weights.items():
            if key == 'plddt':
                wscore += w * r.get('plddt', 0)
            elif key == 'iptm':
                wscore += w * (r.get('iptm', 0) * 100)
            elif key == 'helicity':
                wscore += w * (helicity * 100)
            elif key == 'ptm':
                wscore += w * (r.get('ptm', 0) * 100)
        r['weighted_score'] = wscore / norm if norm else 0

    # Rank
    all_results.sort(key=lambda r: r['weighted_score'], reverse=True)
    top = all_results[:args.top_k]

    # Write final results
    final_dir = os.path.join(output_dir, 'final_results')
    top_k_dir = os.path.join(final_dir, 'topk_structures')
    os.makedirs(top_k_dir, exist_ok=True)

    fieldnames = [
        'rank', 'design_name', 'source_tag',
        'best_model', 'plddt', 'ptm', 'iptm', 'helicity',
        'weighted_score', 'mpnn_score', 'n_mutations', 'design_seq',
        'structure_file',
    ]

    csv_path = os.path.join(final_dir, 'final_topk.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for i, r in enumerate(top, 1):
            tag = r['design_name']
            bm = r.get('best_model', '')

            out_name = f"rank{i:03d}_{tag}.pdb"
            out_path = os.path.join(top_k_dir, out_name)
            if bm != '' and not os.path.exists(out_path):
                cif_path = os.path.join(chai1_out, tag, f"pred.model_idx_{bm}.cif")
                if os.path.exists(cif_path):
                    _cif_to_pdb(cif_path, out_path)

            row = {
                'rank': i,
                'design_name': tag,
                'source_tag': r.get('source_tag', ''),
                'best_model': bm,
                'plddt': f"{r['plddt']:.3f}" if r.get('plddt') else '',
                'ptm': f"{r['ptm']:.4f}" if r.get('ptm') else '',
                'iptm': f"{r['iptm']:.4f}" if r.get('iptm') else '',
                'helicity': f"{r['helicity']:.4f}",
                'weighted_score': f"{r['weighted_score']:.3f}",
                'mpnn_score': f"{r['mpnn_score']:.4f}" if r.get('mpnn_score') is not None else '',
                'n_mutations': r.get('n_mutations', ''),
                'design_seq': r.get('design_seq', ''),
                'structure_file': out_name if os.path.exists(out_path) else '',
            }
            writer.writerow(row)

    print(f'\n[Step 4] Final top {len(top)} MPNN+Chai-1 designs (weighted by {args.weights}):')
    for i, r in enumerate(top, 1):
        print(f'  {i}. {r["design_name"]:30s}  '
              f'pLDDT={r["plddt"]:.3f}  iPTM={r["iptm"]:.4f}  '
              f'helicity={r["helicity"]:.3f}  '
              f'score={r["weighted_score"]:.3f}')
    print(f'\n  Structures: {top_k_dir}/')
    print(f'  CSV: {csv_path}')

    return top


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='RFDiffusion3 Design -> Chai-1 -> ProteinMPNN -> Chai-1',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    # Input / design
    p.add_argument('--pdb', required=True, help='Input PDB file')
    p.add_argument('--design-regions', required=True,
                   help='Design regions: <chain><start>-<end>:<len_min>-<len_max>[;...]')
    p.add_argument('--output', '-o', default='./protocol_out')

    # Step 1: RFD3
    p.add_argument('--num-rounds', type=int, default=1,
                   help='Number of independent RFD3 rounds (default: 1)')
    p.add_argument('--num-designs', type=int, default=50,
                   help='RFD3 designs per round (default: 50)')
    p.add_argument('--top-m', type=int, default=10,
                   help='Top unique RFD3 designs for step 2 Chai-1 (default: 10)')

    # Step 2: Chai-1
    p.add_argument('--top-k', type=int, default=3,
                   help='Top designs after step 2 ranking for MPNN (default: 3)')

    # Step 3: MPNN
    p.add_argument('--mpnn-num-seqs', type=int, default=100,
                   help='Sequences per ProteinMPNN run (default: 100)')
    p.add_argument('--mpnn-temperature', type=float, default=0.1,
                   help='ProteinMPNN sampling temperature (default: 0.1)')
    p.add_argument('--mpnn-design-unfixed', action='store_true',
                   help='Design ALL residues in designed chains (not just --design-regions)')
    p.add_argument('--mpnn-chai1-top', type=int, default=0,
                   help='Number of top MPNN designs for step 4 Chai-1 '
                        '(0 = top half, default: 0)')

    # Slurm / execution
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
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    cfg = load_config(args.config)
    os.makedirs(args.output, exist_ok=True)

    design_regions = parse_design_regions(args.design_regions)
    designed_chains_set = set(r[0] for r in design_regions)
    print(f'Design regions: {design_regions}')
    print(f'Designed chains: {sorted(designed_chains_set)}')

    # ==================================================================
    # Phase 1: RFDiffusion3 backbone generation (N rounds, parallel)
    # ==================================================================
    round_dirs = []
    rfd3_needed = []
    for r in range(1, args.num_rounds + 1):
        rfd3_out = os.path.join(args.output, 'step1_rfd3', f'round_{r}')
        scores_csv = os.path.join(rfd3_out, 'rfd3_scores.csv')
        if args.skip_existing and os.path.exists(scores_csv):
            print(f'[Step 1] Skipping round {r}/{args.num_rounds} -- already completed')
            round_dirs.append(rfd3_out)
        else:
            rfd3_needed.append((r, rfd3_out))

    if rfd3_needed:
        processes = []
        for r, rfd3_out in rfd3_needed:
            print(f'[Step 1] Submitting RFDiffusion3 round {r}/{args.num_rounds} '
                  f'({args.num_designs} designs)')
            cmd = build_rfd3_cmd(args, rfd3_out, args.config)
            proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
            processes.append((r, rfd3_out, proc))

        print(f'\n  Submitted {len(processes)} round(s) -- waiting for all to complete...\n')
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

    # ==================================================================
    # Collect designs, deduplicate, select top M for Chai-1
    # ==================================================================
    all_designs, chain_meta = collect_designs_from_rounds(round_dirs)
    if not all_designs:
        print('ERROR: no designs collected from any round.')
        sys.exit(1)

    unique = deduplicate_designs(all_designs, chain_meta)
    top_m = unique[:args.top_m]
    print(f'\n[Step 1b] Selected top {len(top_m)}/{len(unique)} unique designs for Chai-1')

    if not top_m:
        print('No designs selected. Exiting.')
        return

    # Write unique FASTA files
    unique_fasta_dir = write_unique_fastas(
        top_m, os.path.join(args.output, 'step1_rfd3', 'unique_designs'))

    # ==================================================================
    # Phase 2: Chai-1 prediction + ranking
    # ==================================================================
    print('\n[Step 2] Chai-1 structure prediction')
    chai1_out = os.path.join(args.output, 'step2_chai1')
    results = run_chai1_batch(top_m, unique_fasta_dir, chai1_out, cfg, args)

    if not results:
        print('No Chai-1 results. Exiting.')
        return

    designed_chains = get_designed_chain_ids(chain_meta)
    python_exe = cfg.get('esmif', {}).get('python_executable')

    # Rank step 2 results and save to final_results/step2/
    step2_top = rank_and_save_step2(
        results, top_m, args.output, chai1_out,
        args.top_k, args.weights,
        designed_chains=designed_chains, python_exe=python_exe,
    )

    # ==================================================================
    # Phase 3: ProteinMPNN design on top-K structures
    # ==================================================================
    print('\n[Step 3] ProteinMPNN sequence design')
    mpnn_designs = run_mpnn_phase(
        step2_top, args.output, cfg, args, design_regions)

    if not mpnn_designs:
        print('No MPNN designs generated. Exiting.')
        return

    # ==================================================================
    # Phase 4: Chai-1 validation of MPNN designs + final ranking
    # ==================================================================
    print('\n[Step 4] Chai-1 validation of MPNN designs')
    run_mpnn_chai1(mpnn_designs, args.output, cfg, args, designed_chains)

    print('\nDone.')


if __name__ == '__main__':
    main()
