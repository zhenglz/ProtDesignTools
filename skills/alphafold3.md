# AlphaFold3 Skill

## Overview
AlphaFold3 is a deep learning model for predicting the structure of proteins, protein complexes, and biomolecular interactions with ligands, DNA, RNA, and ions.

## Installation
AlphaFold3 is installed at: `/data_test/home/guoliangzhu/bioapp/alphafold3`

### Environment Setup
```bash
# Navigate to AlphaFold3 directory
cd /data_test/home/guoliangzhu/bioapp/alphafold3

# Source the environment
source source_af3_env.sh

# Or activate the conda environment
conda activate af3

# Check installation
./bin/af3 --help
```

## Basic Usage

**IMPORTANT**: AlphaFold3 requires GPU acceleration and must be run on GPU nodes.

### Direct Execution (Interactive GPU Node)
```bash
# First, get on a GPU node
srun -p 4090 --gres=gpu:1 --time=4:00:00 --pty bash

# Then run AlphaFold3
cd /data_test/home/guoliangzhu/bioapp/alphafold3
./bin/af3 -f <input.fasta> -o <output_directory>

# Example from provided usage
./bin/af3 -f 701.fasta -o 701

# Using the run.sh script (handles environment setup)
./run.sh <input.fasta> <output_directory>
```

### SLURM Submission (Recommended)
```bash
# Create a SLURM script for AlphaFold3
cat > run_af3.slurm << 'EOF'
#!/bin/bash
#SBATCH -J af3_prediction
#SBATCH -p 4090                    # Use 4090 partition for RTX 4090 GPUs
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4GB
#SBATCH --gres=gpu:1               # Request 1 GPU
#SBATCH --time=24:00:00
#SBATCH -o %j.out
#SBATCH -e %j.err

cd /data_test/home/guoliangzhu/bioapp/alphafold3
./bin/af3 -f protein.fasta -o output_AF3
EOF

# Submit the job
sbatch run_af3.slurm
```

### Using the Submission Script
```bash
# Use the provided SLURM submission script
cd /data_test/home/lzzheng/bin
./submit_slurm_gpu.sh "cd /data_test/home/guoliangzhu/bioapp/alphafold3 && ./bin/af3 -f protein.fasta -o output_AF3" 8 4090
```

### Key Parameters
- `-f`: Input FASTA file path
- `-o`: Output directory for prediction results
- `--num_samples`: Number of samples to generate (default: 5)
- `--seed`: Random seed for reproducibility
- `--model_preset`: Model preset (monomer, multimer, etc.)
- `--max_recycles`: Maximum number of recycling iterations
- `--tol`: Tolerance for early stopping

### GPU Requirements
AlphaFold3 requires GPU acceleration. Recommended partitions:
1. **4090 partition**: RTX 4090 GPUs (24GB VRAM) - Best for most predictions
2. **gpu_part**: General GPU partition - May have mixed GPU types
3. **p100 partition**: P100 GPUs (16GB VRAM) - For smaller predictions

**Memory considerations**: Large complexes may require RTX 4090 (24GB VRAM) or multiple GPUs.

## Input Requirements

### FASTA Format
AlphaFold3 automatically recognizes different molecule types from FASTA headers and sequences:

```bash
# Example FASTA file with protein and ligand
cat > complex.fasta << 'EOF'
>701
MANPYERGPNPTDALLEARSGPFSVSEENVSRLSASGFGGGTIYYPRENNTYGAVAISPGYTGTEASIAWLGKRIASHGFVVITIDTITTLDQPDSRAEQLNAALNHMINRASSTVRSRIDSSRLAVMGHSMGGGGSLRLASQRPDLKAAIPLTPWHLNKNWSSVRVPTLIIGADLDTIAPVLTHARPFYNSLPTSISKAYLELDGATHFAPNIPNKIIGKYSVAWLKRFVDNDTRYTQFLCPGPRDGLFGEVEEYRSTCPF
>Substrate 
O=C(O)c2ccc(C(=O)OCCOC(=O)c1ccc(C(=O)OCCO)cc1)cc2
EOF
```

### Supported Molecule Types
AlphaFold3 automatically recognizes:
- **Proteins**: Standard amino acid sequences
- **Ligands**: SMILES strings or small molecule representations
- **DNA**: DNA sequences (ATCG)
- **RNA**: RNA sequences (AUCG)
- **Ions**: Ion specifications

### Input Preparation
```bash
# Create FASTA file with multiple components
cat > protein_ligand.fasta << 'EOF'
>protein_chain_A
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLA
>ligand_1
CC(=O)OC1=CC=CC=C1C(=O)O
>dna_chain_B
ATCGATCGATCG
EOF

# Check sequence composition
python -c "
from Bio import SeqIO
for record in SeqIO.parse('complex.fasta', 'fasta'):
    seq = str(record.seq)
    print(f'{record.id}: {len(seq)} characters')
    print(f'First 50 chars: {seq[:50]}...')
"
```

## Output Structure

### Generated Files
```
output_directory/
├── out/
│   ├── out_model.cif                    # Predicted structure (mmCIF format)
│   ├── out_confidences.json             # Detailed confidence scores
│   ├── out_summary_confidences.json     # Summary confidence scores
│   ├── out_data.json                    # Raw prediction data
│   ├── ranking_scores.csv               # Ranking scores for samples
│   ├── seed-1_sample-0/                 # Sample 0 output
│   ├── seed-1_sample-1/                 # Sample 1 output
│   ├── seed-1_sample-2/                 # Sample 2 output
│   ├── seed-1_sample-3/                 # Sample 3 output
│   └── seed-1_sample-4/                 # Sample 4 output
└── TERMS_OF_USE.md                      # Usage terms
```

### Output Files Description

1. **out_model.cif**: Main predicted structure in mmCIF format
2. **out_summary_confidences.json**: Summary confidence metrics (pTM, iPTM, etc.)
3. **out_confidences.json**: Detailed per-residue confidence scores
4. **ranking_scores.csv**: Scores for ranking different samples
5. **seed-*_sample-*/**: Individual sample outputs

## Output Interpretation

### Confidence Metrics
The key confidence metrics are in `out_summary_confidences.json`:

```json
{
 "chain_iptm": [0.96, 0.8, 0.79],           // iPTM per chain
 "chain_pair_iptm": [[...], [...], [...]],  // iPTM between chain pairs
 "chain_pair_pae_min": [[...], [...]],      // Minimum PAE between chains
 "chain_ptm": [0.92, 0.66, 0.04],           // pTM per chain
 "fraction_disordered": 0.03,               // Fraction disordered residues
 "has_clash": 0.0,                          // Clash indicator
 "iptm": 0.96,                              // Interface pTM (complexes)
 "ptm": 0.93,                               // Predicted TM-score
 "ranking_score": 0.97                      // Overall ranking score
}
```

### Analyzing Confidence Scores
```bash
# Extract key confidence metrics
python -c "
import json

with open('output_AF3_complex/out/out_summary_confidences.json', 'r') as f:
    data = json.load(f)

print('=== AlphaFold3 Confidence Scores ===')
print(f'Overall pTM: {data.get(\"ptm\", 0):.3f}')
print(f'Interface pTM (iPTM): {data.get(\"iptm\", 0):.3f}')
print(f'Ranking score: {data.get(\"ranking_score\", 0):.3f}')
print(f'Fraction disordered: {data.get(\"fraction_disordered\", 0):.3f}')
print(f'Has clash: {data.get(\"has_clash\", 0)}')

# Chain-specific scores
if 'chain_ptm' in data:
    print('\\nChain pTM scores:')
    for i, ptm in enumerate(data['chain_ptm']):
        print(f'  Chain {i}: {ptm:.3f}')

if 'chain_iptm' in data:
    print('\\nChain iPTM scores:')
    for i, iptm in enumerate(data['chain_iptm']):
        print(f'  Chain {i}: {iptm:.3f}')
"

# Check ranking scores
cat output_AF3_complex/out/ranking_scores.csv
```

### Convert CIF to PDB
```bash
# Using the provided conversion script
python /data_test/home/lzzheng/bin/cif2pdb.py out_model.cif output.pdb

# Or using BioPython
python -c "
from Bio.PDB import MMCIFParser, PDBIO
parser = MMCIFParser()
structure = parser.get_structure('model', 'output_AF3_complex/out/out_model.cif')
io = PDBIO()
io.set_structure(structure)
io.save('output_AF3_complex/out/model.pdb')
print('Converted CIF to PDB')
"
```

## Advanced Usage

### Complex Prediction with Multiple Chains
```bash
# FASTA with protein-protein complex
cat > protein_complex.fasta << 'EOF'
>receptor
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLA
>ligand
MKLKKLRLPSLKVLLFIFSLLLACSSPSAQTSPAVQTSPAVQTSPAVQTSPAVQTSPAVQTSP
EOF

# Run prediction
cd /data_test/home/guoliangzhu/bioapp/alphafold3
./bin/af3 -f protein_complex.fasta -o protein_complex_AF3
```

### Protein-Ligand Complex
```bash
# FASTA with protein and small molecule
cat > protein_ligand.fasta << 'EOF'
>enzyme
MANPYERGPNPTDALLEARSGPFSVSEENVSRLSASGFGGGTIYYPRENNTYGAVAISPGYTGTEASIAWLGKRIASHGFVVITIDTITTLDQPDSRAEQLNAALNHMINRASSTVRSRIDSSRLAVMGHSMGGGGSLRLASQRPDLKAAIPLTPWHLNKNWSSVRVPTLIIGADLDTIAPVLTHARPFYNSLPTSISKAYLELDGATHFAPNIPNKIIGKYSVAWLKRFVDNDTRYTQFLCPGPRDGLFGEVEEYRSTCPF
>inhibitor
O=C(O)c2ccc(C(=O)OCCOC(=O)c1ccc(C(=O)OCCO)cc1)cc2
EOF

./bin/af3 -f protein_ligand.fasta -o protein_ligand_AF3
```

### DNA-Protein Complex
```bash
# FASTA with DNA and protein
cat > dna_protein.fasta << 'EOF'
>transcription_factor
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLA
>dna_binding_site
ATCGATCGATCGATCGATCG
EOF

./bin/af3 -f dna_protein.fasta -o dna_protein_AF3
```

## Quality Assessment

### Confidence Score Interpretation
1. **pTM (0-1)**: Predicted TM-score
   - >0.7: Very good prediction
   - 0.5-0.7: Good prediction
   - <0.5: May have issues

2. **iPTM (0-1)**: Interface predicted TM-score (for complexes)
   - >0.5: Good interface prediction
   - <0.5: Poor interface prediction

3. **Ranking score (0-1)**: Overall model quality
   - >0.8: High confidence
   - 0.6-0.8: Medium confidence
   - <0.6: Low confidence

4. **Fraction disordered (0-1)**: Proportion of disordered residues
   - <0.1: Well-structured
   - 0.1-0.3: Some disorder
   - >0.3: Highly disordered

### Quality Check Script
```bash
#!/bin/bash
# quality_check_af3.sh
output_dir=$1

python -c "
import json
import os

conf_file = os.path.join('$output_dir', 'out/out_summary_confidences.json')
if not os.path.exists(conf_file):
    print('Error: Confidence file not found')
    exit(1)

with open(conf_file, 'r') as f:
    data = json.load(f)

ptm = data.get('ptm', 0)
iptm = data.get('iptm', 0)
ranking = data.get('ranking_score', 0)
disordered = data.get('fraction_disordered', 0)
has_clash = data.get('has_clash', 0)

print('=== AlphaFold3 Quality Assessment ===')
print(f'pTM: {ptm:.3f} ({'GOOD' if ptm > 0.7 else 'FAIR' if ptm > 0.5 else 'POOR'})')
print(f'iPTM: {iptm:.3f} ({'GOOD' if iptm > 0.5 else 'POOR'})')
print(f'Ranking score: {ranking:.3f} ({'HIGH' if ranking > 0.8 else 'MEDIUM' if ranking > 0.6 else 'LOW'})')
print(f'Disordered fraction: {disordered:.3f}')
print(f'Has clash: {has_clash}')

# Overall assessment
if ptm > 0.7 and ranking > 0.8:
    print('\\nOverall: HIGH QUALITY PREDICTION')
elif ptm > 0.5 and ranking > 0.6:
    print('\\nOverall: MEDIUM QUALITY PREDICTION')
else:
    print('\\nOverall: LOW QUALITY PREDICTION - consider alternative methods')
"
```

## Integration with Other Tools

### Compare with Experimental Structure
```bash
# Convert AF3 output to PDB
python /data_test/home/lzzheng/bin/cif2pdb.py output_AF3_complex/out/out_model.cif af3_model.pdb

# Align with experimental structure
TM-align experimental.pdb af3_model.pdb -o alignment.txt

# Calculate RMSD
python -c "
from Bio.PDB import PDBParser, Superimposer
parser = PDBParser()
ref = parser.get_structure('ref', 'experimental.pdb')
pred = parser.get_structure('pred', 'af3_model.pdb')

# Superimpose on CA atoms
ref_atoms = [atom for atom in ref.get_atoms() if atom.get_name() == 'CA']
pred_atoms = [atom for atom in pred.get_atoms() if atom.get_name() == 'CA']

min_len = min(len(ref_atoms), len(pred_atoms))
si = Superimposer()
si.set_atoms(ref_atoms[:min_len], pred_atoms[:min_len])
si.apply(pred.get_atoms())

rmsd = si.rms
print(f'RMSD after superposition: {rmsd:.3f} Å')
"
```

### Combine with Chai-1 Predictions
```bash
# Run both AF3 and Chai-1 for comparison
cd /data_test/home/guoliangzhu/bioapp/alphafold3
./bin/af3 -f protein.fasta -o output_AF3

cd /data_test/share/pub_tools/chai-lab
./run.sh protein.fasta output_chai1

# Compare confidence scores
python compare_predictions.py --af3 output_AF3 --chai1 output_chai1
```

## Troubleshooting

### Common GPU Issues
1. **CUDA out of memory**: 
   - Use RTX 4090 partition (24GB VRAM)
   - Reduce sequence length or complexity
   - Use `--num_samples 1` for fewer predictions
   - Check the `run.sh` script for memory optimization flags

2. **No GPU available**:
   - Ensure you're on a GPU node: `srun -p 4090 --gres=gpu:1 --pty bash`
   - Check partition availability: `sinfo -p 4090,gpu_part,p100`
   - Submit via SLURM instead of interactive run

3. **CUDA version mismatch**:
   - AlphaFold3 requires specific CUDA versions
   - Load correct modules: `module load compiler/cuda/12.4` for 4090

4. **FASTA parsing errors**: Ensure proper FASTA format
5. **Output directory exists**: Remove or rename existing output directory

### Debug Commands
```bash
# Test GPU availability
nvidia-smi
echo $CUDA_VISIBLE_DEVICES

# Test AF3 installation on GPU node
srun -p 4090 --gres=gpu:1 --time=1:00:00 --pty bash
cd /data_test/home/guoliangzhu/bioapp/alphafold3
./bin/af3 --help

# Check environment
source source_af3_env.sh
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"

# Test with small sequence
cat > test.fasta << 'EOF'
>test_protein
MKTVRQERLK
EOF
./bin/af3 -f test.fasta -o test_output
```

## Performance Tips

### GPU Resource Optimization
```bash
# For small proteins (<200 residues) - can use p100 partition
./submit_slurm_gpu.sh "cd /data_test/home/guoliangzhu/bioapp/alphafold3 && ./bin/af3 -f small.fasta -o output" 4 p100

# For medium proteins (200-500 residues) - use gpu_part
./submit_slurm_gpu.sh "cd /data_test/home/guoliangzhu/bioapp/alphafold3 && ./bin/af3 -f medium.fasta -o output" 8 gpu_part

# For large complexes (>500 residues or multiple chains) - use 4090 partition
./submit_slurm_gpu.sh "cd /data_test/home/guoliangzhu/bioapp/alphafold3 && ./bin/af3 -f complex.fasta -o output" 8 4090

# Reduce memory usage
./bin/af3 -f complex.fasta -o output --num_samples 1  # Fewer samples
```

### Batch Processing with SLURM Arrays
```bash
#!/bin/bash
#SBATCH -J af3_batch
#SBATCH -p 4090
#SBATCH -a 1-10                    # 10 array jobs
#SBATCH -n 8
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4GB
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00

cd /data_test/home/guoliangzhu/bioapp/alphafold3
fasta_file="sequences/protein_${SLURM_ARRAY_TASK_ID}.fasta"
output_dir="af3_output/run_${SLURM_ARRAY_TASK_ID}"

./bin/af3 -f "$fasta_file" -o "$output_dir"
```

### Memory Management for Large Complexes
The `run.sh` script includes memory optimization flags:
```bash
# Key environment variables for memory management
export XLA_PYTHON_CLIENT_PREALLOCATE=false  # Disable preallocation
export TF_FORCE_UNIFIED_MEMORY=true         # Unified CPU/GPU memory
export XLA_CLIENT_MEM_FRACTION=3.2          # Memory fraction

# Run with memory optimization
cd /data_test/home/guoliangzhu/bioapp/alphafold3
./run.sh large_complex.fasta output_af3
```

## References
- AlphaFold3 installation: `/data_test/home/guoliangzhu/bioapp/alphafold3`
- Main script: `./bin/af3`
- Environment: `source_af3_env.sh`
- CIF to PDB conversion: `/data_test/home/lzzheng/bin/cif2pdb.py`
- Confidence metrics: `out_summary_confidences.json`