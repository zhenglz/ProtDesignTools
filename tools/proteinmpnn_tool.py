#!/usr/bin/env python3
"""
ProteinMPNN Tool - A unified interface for scoring and designing protein sequences
"""

import os
import sys
import argparse
import subprocess as sp
import shutil
import tempfile
import json
import pandas as pd
try:
    import mdtraj as mt
    MDTRAJ_AVAILABLE = True
except ImportError:
    MDTRAJ_AVAILABLE = False
from pathlib import Path


# Configuration
PACKAGE_DPATH = "/data_test/home/lzzheng/apps/ProteinMPNN/ProteinMPNN"
#PYTHON_EXE = "/data_test/home/lzzheng/.conda/envs/sfct/bin/python3.6"
SCORE_PYTHON = f"{PACKAGE_DPATH}/../python_env/proteinmpnn/bin/python"
PYTHON_EXE = f"{PACKAGE_DPATH}/../python_env/proteinmpnn/bin/python"
SCRIPT_DPATH = f"{PACKAGE_DPATH}/../latest/ProteinMPNN"


def argument():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="ProteinMPNN Tool - Score and design protein sequences"
    )

    # Required arguments
    parser.add_argument('-f', '--pdb', required=True, help='Input PDB file path')
    parser.add_argument('-o', '--output', required=True, help='Output directory path')

    # Mode selection
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('--score', action='store_true', help='Score mode')
    mode_group.add_argument('--design', action='store_true', help='Design mode')

    # Scoring options
    parser.add_argument('--seq', help='Sequence to score (FASTA format or sequence string)')
    parser.add_argument('--mut', help='Mutations to score (e.g., "A12G,A13V")')

    # Design options
    parser.add_argument('--positions', help='Positions to design in ${ResidueSeq}${chainid} format (e.g., "21A-25A,145B-150B" or "17A,19A,132A,134B")')
    parser.add_argument('--exclude', help='Positions to exclude from design (keep fixed) in same format as --positions')
    parser.add_argument('-n', type=int, default=100, help='Number of sequences to design (default: 100)')
    parser.add_argument('-t', '--temperature', type=float, default=0.1, help='Sampling temperature (default: 0.1)')

    # Optional
    parser.add_argument('--chain', default='A', help='Chain ID for scoring (default: A)')

    args = parser.parse_args()

    # Validate arguments based on mode
    if args.score:
        if not args.seq and not args.mut:
            print("[INFO] No sequence or mutations provided for scoring. Will score the sequence from PDB.")
    elif args.design:
        if not args.positions:
            parser.error("--design mode requires --positions argument")

    return args


def simple_pdb2fasta(pdb_fpath):
    """Simple PDB to FASTA parser without mdtraj"""
    aa_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    chains = {}
    with open(pdb_fpath) as f:
        for line in f:
            if line.startswith('ATOM'):
                res_num = int(line[22:26].strip())
                chain_id = line[21]
                res_name = line[17:20].strip()

                if chain_id not in chains:
                    chains[chain_id] = {}

                if res_num not in chains[chain_id]:
                    aa_code = aa_dict.get(res_name, 'X')
                    chains[chain_id][res_num] = aa_code

    # Convert to sequences
    seqs = []
    for chain_id in sorted(chains.keys()):
        residues = sorted(chains[chain_id].items())
        seq = ''.join([aa for _, aa in residues])
        seqs.append(seq)

    return "/".join(seqs)


def pdb2fasta(pdb_fpath):
    """Extract sequence from PDB file"""
    if not MDTRAJ_AVAILABLE:
        print("[WARNING] mdtraj not available. Using simple PDB parser.")
        return simple_pdb2fasta(pdb_fpath)

    try:
        pdb = mt.load_pdb(pdb_fpath)
        seq = "/".join(pdb.topology.to_fasta())
        return seq
    except Exception as e:
        print(f"[ERROR] Failed to extract sequence from PDB with mdtraj: {e}")
        print("[INFO] Falling back to simple parser")
        return simple_pdb2fasta(pdb_fpath)


def parse_position_string(position_str):
    """
    Parse position string in ${ResidueSeq}${chainid} format
    Examples: "21A-25A,145B-150B" or "17A,19A,132A,134B"
    Returns list of (residue_index, chain_id) tuples
    """
    positions = []

    if not position_str:
        return positions

    if "-" in position_str:
        # Handle ranges like "21A-25A,145B-150B"
        for segment in position_str.split(","):
            if "-" in segment:
                # Extract range
                start_part = segment.split("-")[0]
                end_part = segment.split("-")[1]

                # Parse residue number and chain
                # Find where number ends and chain starts
                start_num = ''
                start_chain = ''
                for char in start_part:
                    if char.isdigit():
                        start_num += char
                    else:
                        start_chain += char

                end_num = ''
                end_chain = ''
                for char in end_part:
                    if char.isdigit():
                        end_num += char
                    else:
                        end_chain += char

                if start_chain != end_chain:
                    print(f"[WARNING] Chain mismatch in range {segment}: {start_chain} != {end_chain}")
                    continue

                chain = start_chain
                start_idx = int(start_num)
                end_idx = int(end_num)

                for idx in range(start_idx, end_idx + 1):
                    positions.append((idx, chain))
            else:
                # Single residue like "17A"
                # Parse residue number and chain
                num = ''
                chain = ''
                for char in segment:
                    if char.isdigit():
                        num += char
                    else:
                        chain += char

                if num and chain:
                    positions.append((int(num), chain))
    else:
        # Handle single residues like "17A,19A,132A,134B"
        for segment in position_str.split(","):
            # Parse residue number and chain
            num = ''
            chain = ''
            for char in segment:
                if char.isdigit():
                    num += char
                else:
                    chain += char

            if num and chain:
                positions.append((int(num), chain))

    return positions


def parse_residues(design_positions_str, exclude_positions_str, pdb_fpath):
    """
    Parse design and exclusion position strings
    Returns dicts of {chain: [residue_indices]} for design and fixed residues
    """
    # Get all residues from PDB
    all_res = []
    with open(pdb_fpath) as lines:
        for x in lines:
            if x.startswith('ATOM'):
                res_idx = int(x[22:26].strip())
                chain_id = x[21]
                res = (res_idx, chain_id)
                if res not in all_res:
                    all_res.append(res)

    # Parse design positions
    design_res = parse_position_string(design_positions_str)

    # Parse exclusion positions
    exclude_res = parse_position_string(exclude_positions_str)

    # Remove excluded residues from design list
    design_res = [res for res in design_res if res not in exclude_res]

    print(f"[INFO] Selected residues for design: {design_res}")
    if exclude_res:
        print(f"[INFO] Excluded residues (kept fixed): {exclude_res}")

    # Calculate fixed residues (all residues not in design_res, plus excluded residues)
    fix_res = [x for x in all_res if x not in design_res]

    # Convert to dict format
    design_res_dict = {}
    for (idx, chain) in design_res:
        if chain not in design_res_dict:
            design_res_dict[chain] = [idx]
        else:
            design_res_dict[chain].append(idx)

    fix_res_dict = {}
    for (idx, chain) in fix_res:
        if chain not in fix_res_dict:
            fix_res_dict[chain] = [idx]
        else:
            fix_res_dict[chain].append(idx)

    return design_res_dict, fix_res_dict


def create_fasta_file(seq_str, output_dir):
    """Create a temporary FASTA file from sequence string"""
    fasta_path = os.path.join(output_dir, "input_seq.fasta")

    # Check if seq_str is a file path
    if os.path.exists(seq_str):
        shutil.copy(seq_str, fasta_path)
        return fasta_path

    # Otherwise treat as sequence string
    with open(fasta_path, 'w') as f:
        f.write(">input_sequence\n")
        f.write(seq_str.replace("/", "\n").replace("\\", "\n"))

    return fasta_path


def apply_mutations(wt_seq, mut_str):
    """
    Apply mutations to wild-type sequence
    mut_str format: "A12G,A13V" where A=original, 12=position, G=mutation
    """
    # Split sequence if it contains chain separators
    if "/" in wt_seq:
        chains = wt_seq.split("/")
    else:
        chains = [wt_seq]

    # Parse mutations
    mutations = mut_str.split(",")

    # Apply mutations to each chain
    mutated_chains = []
    for chain_idx, chain_seq in enumerate(chains):
        chain_id = chr(65 + chain_idx)  # A, B, C, ...
        chain_mut_seq = list(chain_seq)

        for mut in mutations:
            if len(mut) < 3:
                continue

            orig_aa = mut[0]
            pos = int(mut[1:-1])
            new_aa = mut[-1]

            # Check if mutation applies to this chain
            # Simple heuristic: if position is within chain length
            if 1 <= pos <= len(chain_seq):
                # Also check original AA matches
                if chain_seq[pos-1] == orig_aa:
                    chain_mut_seq[pos-1] = new_aa

        mutated_chains.append(''.join(chain_mut_seq))

    return "/".join(mutated_chains)


def run_design(pdb_fpath, positions, exclude_positions, output_dir, num_seqs=100, temperature=0.1):
    """Run ProteinMPNN design mode"""
    os.makedirs(output_dir, exist_ok=True)

    # Create input directory
    input_dir = os.path.join(output_dir, 'input')
    os.makedirs(input_dir, exist_ok=True)
    shutil.copy(pdb_fpath, os.path.join(input_dir, 'protein.pdb'))

    # Parse positions
    design_res_dict, fix_res_dict = parse_residues(positions, exclude_positions, pdb_fpath)

    # Prepare fixed positions list
    fix_res_list = []
    for chain in fix_res_dict.keys():
        fix_res_list.append(" ".join([str(x) for x in fix_res_dict[chain]]))

    # Create run script
    script_path = os.path.join(output_dir, 'run_design.sh')

    with open(script_path, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write(f'''
folder_with_pdbs={input_dir}
output_dir={output_dir}

path_for_parsed_chains=$output_dir/parsed_pdbs.jsonl
path_for_assigned_chains=$output_dir/assigned_pdbs.jsonl
path_for_fixed_positions=$output_dir/fixed_pdbs.jsonl

chains_to_design="{" ".join(list(fix_res_dict.keys()))}"

# define fixed residues
fixed_positions="{",".join(fix_res_list)}"

{PYTHON_EXE} {PACKAGE_DPATH}/vanilla_proteinmpnn/helper_scripts/parse_multiple_chains.py \\
             --input_path=$folder_with_pdbs \\
             --output_path=$path_for_parsed_chains

{PYTHON_EXE} {PACKAGE_DPATH}/vanilla_proteinmpnn/helper_scripts/assign_fixed_chains.py \\
             --input_path=$path_for_parsed_chains \\
             --output_path=$path_for_assigned_chains \\
             --chain_list "$chains_to_design"

{PYTHON_EXE} {PACKAGE_DPATH}/vanilla_proteinmpnn/helper_scripts/make_fixed_positions_dict.py \\
             --input_path=$path_for_parsed_chains \\
             --output_path=$path_for_fixed_positions \\
             --chain_list "$chains_to_design" \\
             --position_list "$fixed_positions"

# run design
{PYTHON_EXE} {PACKAGE_DPATH}/vanilla_proteinmpnn/protein_mpnn_run.py \\
        --jsonl_path $path_for_parsed_chains \\
        --chain_id_jsonl $path_for_assigned_chains \\
        --fixed_positions_jsonl $path_for_fixed_positions \\
        --out_folder $output_dir \\
        --num_seq_per_target {num_seqs} \\
        --sampling_temp {temperature} \\
        --batch_size 1
''')

    # Make script executable and run
    os.chmod(script_path, 0o755)

    print(f"[INFO] Running design script: {script_path}")
    result = sp.run(f"bash {script_path}", shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"[ERROR] Design failed: {result.stderr}")
        return False

    return True


def run_score(pdb_fpath, seq_str, output_dir, chain_id='A'):
    """Run ProteinMPNN score mode"""
    os.makedirs(output_dir, exist_ok=True)

    # Create FASTA file
    fasta_path = create_fasta_file(seq_str, output_dir)

    # Run scoring
    cmd = [
        SCORE_PYTHON, f"{SCRIPT_DPATH}/protein_mpnn_run.py",
        "--path_to_fasta", fasta_path,
        "--pdb_path", pdb_fpath,
        "--pdb_path_chains", chain_id,
        "--out_folder", output_dir,
        "--num_seq_per_target", "5",
        "--sampling_temp", "0.1",
        "--score_only", "1",
        "--seed", "13",
        "--batch_size", "1"
    ]

    print(f"[INFO] Running scoring command: {' '.join(cmd)}")
    result = sp.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"[ERROR] Scoring failed: {result.stderr}")
        return None

    return output_dir


def parse_design_output(output_dir, wt_seq):
    """Parse design output and return DataFrame"""
    seqs_dir = os.path.join(output_dir, 'seqs')

    if not os.path.exists(seqs_dir):
        print(f"[ERROR] No seqs directory found in {output_dir}")
        return None

    # Find the output file
    seq_files = [f for f in os.listdir(seqs_dir) if f.endswith('.fa') or f.endswith('.fasta')]
    if not seq_files:
        print(f"[ERROR] No sequence files found in {seqs_dir}")
        return None

    seq_file = os.path.join(seqs_dir, seq_files[0])

    # Parse the file
    data = []
    with open(seq_file) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            header = lines[i].strip()
            i += 1
            if i < len(lines):
                seq = lines[i].strip()
                i += 1

                # Parse score from header
                if 'sample=' not in header:
                    name = 'wild_type'
                    score = float(header.split(",")[1].split("=")[1])
                    recovery = 1.0
                else:
                    # For designed sequences
                    name = f"design_{len(data)}"
                    score = float(header.split(",")[2].split("=")[1])
                    recovery = float(header.split(",")[3].split("=")[1])

                data.append([name, score, recovery, seq])
        else:
            i += 1

    df = pd.DataFrame(data, columns=['name', 'score', 'recovery', 'sequence'])

    # Separate wild-type from designed sequences
    wt_mask = df['name'] == 'wild_type'
    wt_df = df[wt_mask]
    design_df = df[~wt_mask]

    # Sort designed sequences in ascending order (smaller scores first)
    design_df = design_df.sort_values('score', ascending=True, ignore_index=True)

    # Reset design sequence numbering
    for i in range(len(design_df)):
        design_df.at[i, 'name'] = f"design_{i}"

    # Combine with wild-type at the top
    df = pd.concat([wt_df, design_df], ignore_index=True)

    # Save to CSV
    csv_path = os.path.join(output_dir, 'design_results.csv')
    df.to_csv(csv_path, index=False)

    return df


def parse_score_output(output_dir):
    """Parse score output and return score"""
    seqs_dir = os.path.join(output_dir, 'seqs')

    if not os.path.exists(seqs_dir):
        print(f"[ERROR] No seqs directory found in {output_dir}")
        return None

    # Find the output file
    seq_files = [f for f in os.listdir(seqs_dir) if f.endswith('.fa') or f.endswith('.fasta')]
    if not seq_files:
        print(f"[ERROR] No sequence files found in {seqs_dir}")
        return None

    seq_file = os.path.join(seqs_dir, seq_files[0])

    # Parse the file - should have just one sequence (the scored one)
    with open(seq_file) as f:
        lines = f.readlines()

    for i in range(len(lines)):
        if lines[i].startswith('>'):
            header = lines[i].strip()
            # Extract score
            if 'score=' in header:
                score = float(header.split('score=')[1].split(',')[0])
                return score

    return None


def main():
    args = argument()

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    if args.score:
        print("[INFO] Running in SCORE mode")

        # Get sequence to score
        if args.seq:
            seq_to_score = args.seq
        elif args.mut:
            # Apply mutations to wild-type sequence
            wt_seq = pdb2fasta(args.pdb)
            if not wt_seq:
                print("[ERROR] Failed to extract sequence from PDB")
                sys.exit(1)
            seq_to_score = apply_mutations(wt_seq, args.mut)
            print(f"[INFO] Applied mutations: {args.mut}")
            print(f"[INFO] Resulting sequence: {seq_to_score}")
        else:
            # Score the sequence from PDB
            seq_to_score = pdb2fasta(args.pdb)
            if not seq_to_score:
                print("[ERROR] Failed to extract sequence from PDB")
                sys.exit(1)
            print(f"[INFO] Scoring sequence from PDB: {seq_to_score}")

        # Run scoring
        score_dir = os.path.join(args.output, "score_results")
        result = run_score(args.pdb, seq_to_score, score_dir, args.chain)

        if result:
            score = parse_score_output(score_dir)
            if score is not None:
                print(f"\n[RESULT] Score: {score}")

                # Save result to file
                result_file = os.path.join(args.output, "score_result.txt")
                with open(result_file, 'w') as f:
                    f.write(f"PDB: {args.pdb}\n")
                    if args.seq:
                        f.write(f"Sequence: {args.seq}\n")
                    if args.mut:
                        f.write(f"Mutations: {args.mut}\n")
                    f.write(f"Scored sequence: {seq_to_score}\n")
                    f.write(f"Score: {score}\n")
                print(f"[INFO] Results saved to {result_file}")
            else:
                print("[ERROR] Failed to parse score from output")
        else:
            print("[ERROR] Scoring failed")

    elif args.design:
        print("[INFO] Running in DESIGN mode")

        # Run design
        design_success = run_design(
            args.pdb,
            args.positions,
            args.exclude,
            args.output,
            num_seqs=args.n,
            temperature=args.temperature
        )

        if design_success:
            # Parse output
            wt_seq = pdb2fasta(args.pdb)
            df = parse_design_output(args.output, wt_seq)

            if df is not None:
                print(f"\n[RESULT] Designed {len(df)-1} sequences (plus wild-type)")
                print("\nTop 10 sequences (wild-type first, then best designs with lowest scores):")
                print(df.head(10).to_string(index=False))

                # Save full results
                result_file = os.path.join(args.output, "design_results.csv")
                df.to_csv(result_file, index=False)
                print(f"\n[INFO] Full results saved to {result_file}")

                # Also save a summary
                summary_file = os.path.join(args.output, "design_summary.txt")
                with open(summary_file, 'w') as f:
                    f.write(f"PDB: {args.pdb}\n")
                    f.write(f"Positions to design: {args.positions}\n")
                    f.write(f"Number of sequences: {args.n}\n")
                    f.write(f"Temperature: {args.temperature}\n")
                    f.write(f"\nTop 10 sequences (wild-type first, then best designs with lowest scores):\n")
                    f.write(df.head(10).to_string(index=False))
                print(f"[INFO] Summary saved to {summary_file}")
            else:
                print("[ERROR] Failed to parse design output")
        else:
            print("[ERROR] Design failed")


if __name__ == "__main__":
    main()
