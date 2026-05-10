#!/usr/bin/env python3
"""
Shared utilities for design protocols: multi-round RFDiffusion3,
sequence deduplication, and unique FASTA output.

Used by design_validate.py and fragment_design_mpnn.py.
"""

import csv, json, os, shutil, sys
from collections import OrderedDict
from pathlib import Path


# ---------------------------------------------------------------------------
# Chain / sequence helpers
# ---------------------------------------------------------------------------

def get_designed_chain_ids(chain_meta):
    """Return the set of chain IDs that have any designed segments."""
    if not chain_meta:
        return set()
    return set(seg['chain'] for seg in chain_meta if seg.get('is_designed'))


def read_fasta_chains(fasta_path):
    """Parse a multi-chain FASTA file.

    Returns OrderedDict mapping chain_id -> sequence, preserving file order.
    Chain ID is extracted from the header: >protein|tag_CHAIN
    """
    chains = OrderedDict()
    cur_id = None
    cur_seq = []

    try:
        with open(fasta_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    # Flush previous chain
                    if cur_id is not None:
                        chains[cur_id] = ''.join(cur_seq)
                    # Extract chain ID: ">protein|design_01_A ..." -> "A"
                    header = line[1:]  # strip '>'
                    parts = header.split('|')
                    if len(parts) >= 2:
                        tag_part = parts[1].split()[0]  # e.g. "design_01_A"
                        # Last underscore segment is chain ID
                        if '_' in tag_part:
                            cur_id = tag_part.rsplit('_', 1)[-1]
                        else:
                            cur_id = tag_part
                    else:
                        cur_id = header.split()[0]
                    cur_seq = []
                else:
                    cur_seq.append(line)
        if cur_id is not None:
            chains[cur_id] = ''.join(cur_seq)
    except Exception as e:
        print(f"  Warning: could not read FASTA {fasta_path}: {e}")
        return OrderedDict()

    return chains


def build_dedup_key(fasta_path, designed_chains):
    """Build a deduplication key from designed-region sequences.

    Reads the multi-chain FASTA, extracts sequences for designed chains only,
    and concatenates them (sorted by chain ID for deterministic ordering).

    Returns a string key, or None if no designed chains were found.
    """
    chains = read_fasta_chains(fasta_path)
    if not chains:
        return None

    if designed_chains:
        # Only designed chains matter for dedup
        parts = []
        for cid in sorted(designed_chains):
            if cid in chains:
                parts.append(f"{cid}:{chains[cid]}")
        if parts:
            return "|".join(parts)
        # All designed chains missing from FASTA — fall through

    # Fallback: use all chains
    parts = [f"{cid}:{seq}" for cid, seq in chains.items()]
    return "|".join(parts) if parts else None


# ---------------------------------------------------------------------------
# CSV score helpers
# ---------------------------------------------------------------------------

_COMPOSITE_FIELDS = [
    "max_ca_deviation", "n_chainbreaks", "n_clashing",
    "radius_of_gyration", "alanine_content", "glycine_content",
]


def composite_rank(row):
    """Compute composite score from a CSV row dict (lower = better).

    Same formula as rfdiffusion3_tool.py.
    """
    return (
        int(float(row.get('n_chainbreaks', 0))) * 100 +
        int(float(row.get('n_clashing', 0))) * 10 +
        abs(float(row.get('radius_of_gyration', 0)) - 12) * 0.1 +
        float(row.get('alanine_content', 0)) * 10 +
        abs(float(row.get('glycine_content', 0)) - 0.05) * 10 +
        float(row.get('max_ca_deviation', 0))
    )


# ---------------------------------------------------------------------------
# Multi-round collection
# ---------------------------------------------------------------------------

def collect_designs_from_rounds(round_dirs):
    """Collect all designs from multiple RFD3 output directories.

    Args:
        round_dirs: list of directory paths (each from one RFD3 round)

    Returns:
        (all_designs, chain_meta) where:
          all_designs: list of dicts with keys:
            design_name, round_dir, fasta_path, score, dedup_key
          chain_meta: the _chain_meta from the first round's rf3.json
            (all rounds share the same meta)
    """
    all_designs = []
    chain_meta = []

    for rdir in round_dirs:
        csv_path = os.path.join(rdir, 'rfd3_scores.csv')
        if not os.path.exists(csv_path):
            print(f"  WARNING: no rfd3_scores.csv in {rdir}, skipping round")
            continue

        fasta_dir = os.path.join(rdir, 'top_designs_fastas')

        # Load chain_meta from the separate metadata file (all rounds share the same)
        if not chain_meta:
            meta_path = os.path.join(rdir, 'rf3_chain_meta.json')
            if os.path.exists(meta_path):
                try:
                    with open(meta_path) as f:
                        chain_meta = json.load(f)
                except Exception:
                    pass

        # Read scores CSV
        with open(csv_path) as f:
            for row in csv.DictReader(f):
                name = row.get('design_name', '')
                fasta_path = os.path.join(fasta_dir, f'{name}.fasta')

                if not os.path.exists(fasta_path):
                    continue

                row['_round_dir'] = rdir
                row['_fasta_path'] = fasta_path
                row['_score'] = composite_rank(row)
                all_designs.append(row)

    return all_designs, chain_meta


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------

def deduplicate_designs(all_designs, chain_meta):
    """Group designs by designed-region sequences, keep best per group.

    Args:
        all_designs: list of design dicts from collect_designs_from_rounds()
        chain_meta: _chain_meta list from rf3.json

    Returns:
        list of unique design dicts, sorted by composite score (best first)
    """
    if not all_designs:
        return []

    designed_chains = get_designed_chain_ids(chain_meta)

    # Group by dedup key
    groups = OrderedDict()
    for d in all_designs:
        key = build_dedup_key(d['_fasta_path'], designed_chains)
        if key is None:
            key = d['design_name']  # unique fallback — keep all
        if key not in groups:
            groups[key] = []
        groups[key].append(d)

    # Keep best (lowest score) per group
    unique = []
    for key, group in groups.items():
        group.sort(key=lambda d: d['_score'])
        best = group[0]
        best['_dedup_group_size'] = len(group)
        unique.append(best)

    # Sort by score (best first)
    unique.sort(key=lambda d: d['_score'])

    n_total = len(all_designs)
    n_unique = len(unique)
    if n_total > n_unique:
        print(f'  Deduplication: {n_total} designs -> {n_unique} unique '
              f'({n_total - n_unique} duplicates removed)')
        if designed_chains:
            print(f'  Designed chains used for dedup: {sorted(designed_chains)}')

    return unique


# ---------------------------------------------------------------------------
# Unique FASTA output
# ---------------------------------------------------------------------------

def write_unique_fastas(unique_designs, output_dir, top_n=0):
    """Write curated FASTA files for unique designs.

    Creates output_dir/unique_fastas/ with one FASTA per unique design,
    and output_dir/unique_designs.fasta with all combined.

    Each design gets a new tag: design_01, design_02, ...
    The tag is stored as d['_tag'] on each design dict.

    Args:
        unique_designs: list of unique design dicts (sorted by score)
        output_dir: directory to write FASTA files into
        top_n: if > 0, keep only top N designs

    Returns:
        fasta_dir path
    """
    os.makedirs(output_dir, exist_ok=True)

    if top_n > 0:
        unique_designs = unique_designs[:top_n]

    fasta_dir = os.path.join(output_dir, 'unique_fastas')
    os.makedirs(fasta_dir, exist_ok=True)

    combined_path = os.path.join(output_dir, 'unique_designs.fasta')
    combined_fh = open(combined_path, 'w')

    for rank, d in enumerate(unique_designs, 1):
        design_tag = f'design_{rank:02d}'
        d['_tag'] = design_tag
        d['_unique_fasta'] = os.path.join(fasta_dir, f'{design_tag}.fasta')
        src_fasta = d['_fasta_path']

        # Read original multi-chain FASTA
        try:
            with open(src_fasta) as f:
                orig_content = f.read()
        except Exception as e:
            print(f"  WARNING: cannot read {src_fasta}: {e}")
            continue

        # Replace design tag in headers to reflect new rank
        # Old header: >protein|design_XX_A ca_dev=...
        # New header: >protein|design_NN_A ca_dev=...
        import re
        new_content = re.sub(
            r'(>[A-Za-z]+\|)design_\d+_',
            rf'\1{design_tag}_',
            orig_content
        )

        # Write per-design FASTA
        with open(d['_unique_fasta'], 'w') as f:
            f.write(new_content)

        # Append to combined FASTA
        combined_fh.write(new_content)

    combined_fh.close()
    print(f'  Unique FASTA: {len(unique_designs)} designs -> {fasta_dir}/')
    print(f'  Combined FASTA: {combined_path}')

    return fasta_dir
