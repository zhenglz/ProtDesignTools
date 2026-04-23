# ProteinMPNN Skill

## Overview
ProteinMPNN is a neural network for protein sequence design that generates amino acid sequences conditioned on a protein backbone structure.

## Installation

### Standard Installation
```bash
# Clone the repository
git clone https://github.com/dauparas/ProteinMPNN.git
cd ProteinMPNN

# Install dependencies
pip install torch
pip install -r requirements.txt
```

### Local Installation (Specific to this system)
ProteinMPNN is installed at: `/data_test/home/lzzheng/apps/ProteinMPNN/ProteinMPNN`

Python environment: `/data_test/home/lzzheng/apps/ProteinMPNN/ProteinMPNN/envs/proteinmpnn`

To activate the environment:
```bash
source /data_test/home/lzzheng/apps/ProteinMPNN/ProteinMPNN/envs/proteinmpnn/bin/activate
```

Or use the Python executable directly:
```bash
/data_test/home/lzzheng/apps/ProteinMPNN/ProteinMPNN/envs/proteinmpnn/bin/python
```

## Basic Usage

### Using Local Installation
```bash
# Using the wrapper script
cd /data_test/home/lzzheng/apps/ProteinMPNN
python run_ProteinMPNN.py -f <input_pdb> -o <output_dir> -n <num_sequences> -t <temperature>

# Example
python run_ProteinMPNN.py -f protein.pdb -o ./output -n 200 -t 0.1

# With specific mutation sites
python run_ProteinMPNN.py -f protein.pdb -o ./output --mut "12A,13A,14A"

# Scoring only (no sequence design)
./run_score_only.sh <input_pdb> <input_fasta> <output_dir> <chain_id>
```

### Direct Usage (Standard ProteinMPNN)
```bash
# Basic usage with a single PDB file
cd /data_test/home/lzzheng/apps/ProteinMPNN/ProteinMPNN/vanilla_proteinmpnn
python protein_mpnn_run.py \
  --pdb_path <input_pdb_file> \
  --out_folder <output_directory> \
  --num_seq_per_target <number_of_sequences> \
  --batch_size <batch_size>

# Example with specific parameters
python protein_mpnn_run.py \
  --pdb_path input.pdb \
  --out_folder ./output \
  --num_seq_per_target 100 \
  --batch_size 10 \
  --sampling_temp "0.1" \
  --seed 37
```

### Key Parameters

#### Local Wrapper Script (`run_ProteinMPNN.py`)
- `-f`: Input PDB file path
- `-o`: Output directory for generated sequences
- `-n`: Number of sequences to generate (default: 200)
- `-t`: Sampling temperature (default: 0.1, lower = more conservative)
- `--mut`: Mutation sites in format "12A,13A,14A" (residue number + chain)

#### Standard ProteinMPNN (`protein_mpnn_run.py`)
- `--pdb_path`: Input PDB file path
- `--out_folder`: Output directory for generated sequences
- `--num_seq_per_target`: Number of sequences to generate per structure
- `--batch_size`: Batch size for inference
- `--sampling_temp`: Sampling temperature (lower = more conservative)
- `--seed`: Random seed for reproducibility
- `--ca_only`: Use only Cα atoms (for low-resolution structures)
- `--chain_id`: Specify which chains to design (comma-separated)

## Input Requirements

### Input Files
1. **PDB file**: Standard protein structure file (.pdb or .cif format)
2. **Optional**: JSON file specifying designable positions

### Input Preparation
```bash
# Ensure PDB file is clean (remove heteroatoms, waters)
grep -E "^ATOM|^HETATM" input.pdb | grep -v "HOH" > clean.pdb

# For local installation
cd /data_test/home/lzzheng/apps/ProteinMPNN
python run_ProteinMPNN.py -f clean.pdb -o ./output

# For multi-chain proteins with standard ProteinMPNN
cd /data_test/home/lzzheng/apps/ProteinMPNN/ProteinMPNN/vanilla_proteinmpnn
python protein_mpnn_run.py --pdb_path complex.pdb --chain_id "A,B"
```

## Output Structure

### Local Wrapper Script Output
```
output_directory/
├── input/
│   └── protein.pdb                    # Copied input PDB
├── run_task.sh                        # Generated bash script
├── seqs/
│   ├── parsed_pdbs.json               # Parsed structure info
│   └── protein_<timestamp>_*.fa       # FASTA sequences
├── protein.fa                         # Combined FASTA file
└── log.txt                            # Run log
```

### Standard ProteinMPNN Output
```
output_directory/
├── seqs/
│   ├── <pdb_basename>_<timestamp>_<model_number>.fa  # FASTA sequences
│   └── parsed_pdbs.json                              # Parsed structure info
├── <pdb_basename>.fa                                 # Combined FASTA file
└── log.txt                                           # Run log
```

### Output Files Description

1. **FASTA files** (`*.fa`):
   - Contains designed sequences in FASTA format
   - Header format: `>design_<index>_score=<score>`
   - Sequence: Amino acid string

2. **run_task.sh**:
   - Generated bash script that runs ProteinMPNN
   - Can be executed manually if needed
   - Contains all parameters and paths

3. **JSON file** (`parsed_pdbs.json`):
   - Contains parsed structure information
   - Includes chain IDs, residue numbers, and amino acid types

4. **Log file** (`log.txt`):
   - Contains runtime information and errors

## Output Interpretation

### Sequence Quality Metrics
1. **Sequence score**: Lower scores indicate better sequences (negative log-likelihood)
   - Typical range: -2.0 to 2.0
   - Lower is better (more native-like)

2. **Sequence diversity**:
   - Compare multiple generated sequences
   - High diversity = more exploration of sequence space

3. **Per-position probabilities**:
   - Available in detailed output mode
   - Shows confidence for each amino acid position

### Example Output Analysis
```bash
# Extract top scoring sequences
grep "score=" output/seqs/*.fa | sort -t'=' -k2 -n | head -10

# Count unique sequences
cat output/*.fa | grep -v "^>" | sort | uniq | wc -l

# Calculate sequence similarity
python -c "
from Bio import pairwise2
from Bio import SeqIO

seqs = [str(rec.seq) for rec in SeqIO.parse('output/designs.fa', 'fasta')]
for i in range(min(5, len(seqs))):
    for j in range(i+1, min(5, len(seqs))):
        score = pairwise2.align.globalxx(seqs[i], seqs[j])[0].score
        similarity = score / max(len(seqs[i]), len(seqs[j]))
        print(f'Seq {i} vs {j}: {similarity:.2%}')
"
```

## Advanced Usage

### Generated run_task.sh Example
The wrapper script generates a `run_task.sh` file. Example content:
```bash
#!/bin/bash

folder_with_pdbs=/path/to/output/input
output_dir=/path/to/output

cd /data_test/home/lzzheng/apps/ProteinMPNN/ProteinMPNN/vanilla_proteinmpnn

PYTHON_EXE=/data_test/home/lzzheng/apps/ProteinMPNN/ProteinMPNN/envs/proteinmpnn/bin/python

$PYTHON_EXE protein_mpnn_run.py \
    --pdb_path "$folder_with_pdbs" \
    --out_folder "$output_dir" \
    --num_seq_per_target 200 \
    --batch_size 10 \
    --sampling_temp "0.1" \
    --seed 37 \
    --suppress_print 0
```

### Designing Specific Regions
```bash
# Create a JSON file for designable positions
cat > designable_positions.json << 'EOF'
{
  "A": {
    "designable": [10, 11, 12, 13, 14],
    "fixed": [1, 2, 3, 4, 5]
  }
}
EOF

# Run with designable positions using standard ProteinMPNN
cd /data_test/home/lzzheng/apps/ProteinMPNN/ProteinMPNN/vanilla_proteinmpnn
python protein_mpnn_run.py \
  --pdb_path input.pdb \
  --out_folder ./output \
  --designable_positions designable_positions.json
```

### Using Different Models
```bash
# Use CA-only model for low-resolution structures
python protein_mpnn_run.py \
  --pdb_path input.pdb \
  --ca_only \
  --model_name "ca_model_weights"

# Use sidechain model for high-resolution structures
python protein_mpnn_run.py \
  --pdb_path input.pdb \
  --model_name "sidechain_model_weights"
```

## Troubleshooting

### Local Installation Issues
1. **Python environment not found**: Ensure the environment exists at `/data_test/home/lzzheng/apps/ProteinMPNN/ProteinMPNN/envs/proteinmpnn`
2. **Permission errors**: Check execute permissions on scripts
3. **Path errors**: Run from `/data_test/home/lzzheng/apps/ProteinMPNN` directory

### Common Issues
1. **PDB parsing errors**: Ensure PDB file follows standard format
2. **Memory issues**: Reduce batch size for large proteins
3. **GPU out of memory**: Use `--batch_size 1` or CPU mode

### Debug Commands
```bash
# Check local installation
cd /data_test/home/lzzheng/apps/ProteinMPNN
ls -la

# Check Python environment
/data_test/home/lzzheng/apps/ProteinMPNN/ProteinMPNN/envs/proteinmpnn/bin/python --version

# Check PDB file validity
python -c "
from Bio.PDB import PDBParser
parser = PDBParser()
structure = parser.get_structure('test', 'input.pdb')
print(f'Chains: {[chain.id for chain in structure.get_chains()]}')
print(f'Residues: {sum(1 for _ in structure.get_residues())}')
"

# Test with small batch size using local wrapper
cd /data_test/home/lzzheng/apps/ProteinMPNN
python run_ProteinMPNN.py -f input.pdb -o test_output -n 5
```

## Integration with Other Tools

### Downstream Analysis
```bash
# Fold sequences with AlphaFold2
for seq_file in output/seqs/*.fa; do
    python alphafold_run.py --fasta_path "$seq_file" --output_dir ./alphafold_results
done

# Score sequences with ESMFold
python esmfold_run.py --fasta_file output/designs.fa --output_dir ./esmfold_results
```

### Scoring Existing Structures with ProteinMPNN
ProteinMPNN can also be used to score existing structures (not just design sequences). This is used in the chimeric design pipeline:

```bash
# Using the scoring script from the chimeric design pipeline
cd /data_test/home/lzzheng/apps/ProtDesignTools
python tools/chimeric_design.py --score-only --fasta input.fasta --output ./results

# Direct scoring using the run_score_only.sh script
cd /data_test/home/lzzheng/apps/ProteinMPNN
./run_score_only.sh <input_pdb> <input_fasta> <output_dir> <chain_id>
```

#### Example for scoring a single structure:
```bash
cd /data_test/home/lzzheng/apps/ProteinMPNN
./run_score_only.sh protein.pdb protein.fasta ./mpnn_scores A
```

#### Output from scoring:
- NPZ files with score arrays
- Best score represents the negative log-likelihood of the sequence given the structure
- Lower scores indicate better compatibility

### Integration with Chimeric Design Pipeline
The chimeric design pipeline (`tools/chimeric_design.py`) uses ProteinMPNN to score chimeric structures:
1. Generates chimeric sequences by swapping loops
2. Predicts structures with Chai-1
3. Scores structures with ProteinMPNN
4. Ranks chimeras by combined scores

```bash
cd /data_test/home/lzzheng/apps/ProtDesignTools
python tools/chimeric_design.py --fasta PH20M3.fasta --msa PH20M3_uniprot.a3m \
  --loops "270-284;68-74" --output ./chimeric_results --submit-slurm
```

## References
- [ProteinMPNN GitHub](https://github.com/dauparas/ProteinMPNN)
- [ProteinMPNN Paper](https://www.science.org/doi/10.1126/science.add2187)