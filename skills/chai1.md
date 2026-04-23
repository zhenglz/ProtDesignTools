# Chai-1 Skill

## Overview
Chai-1 is a protein structure prediction tool similar to AlphaFold2, designed for predicting protein structures from amino acid sequences. It can optionally use MSA (Multiple Sequence Alignment) information for improved predictions.

## Installation
Chai-1 is installed at: `/data_test/share/pub_tools/chai-lab/`

### Environment Setup
```bash
# Navigate to Chai-1 directory
cd /data_test/share/pub_tools/chai-lab

# The tool uses a Python environment at:
# /data_test/share/pub_tools/chai-lab/envs/chai1/bin/python

# Check installation
ls -la run.sh
```

## Basic Usage

### Direct Execution (Interactive)
```bash
# Basic usage
cd /data_test/share/pub_tools/chai-lab
./run.sh <input.fasta> <output_directory> [msa.a3m]

# Example without MSA
./run.sh input.fasta ./output_chai1

# Example with MSA
./run.sh input.fasta ./output_chai1 msa.a3m
```

### SLURM Submission (Recommended for GPU nodes)
```bash
# Using the submit_slurm_gpu.sh script
cd /data_test/home/lzzheng/bin
./submit_slurm_gpu.sh "<command>" <ncpus> <partition>

# Example: Run Chai-1 on 4090 partition with 4 CPUs
./submit_slurm_gpu.sh "/data_test/share/pub_tools/chai-lab/run.sh input.fasta ./output_chai1" 4 4090

# Example: Run on gpu_part partition
./submit_slurm_gpu.sh "/data_test/share/pub_tools/chai-lab/run.sh input.fasta ./output_chai1" 4 gpu_part

# Note: No need to cd to chai-lab directory first - run.sh finds its own path
```

### Key Parameters
- `input.fasta`: Input protein sequence in FASTA format
- `output_directory`: Directory for prediction results
- `msa.a3m` (optional): MSA file in A3M format for improved accuracy
- `--use-msa-server`: Flag to use MSA server (automatically added when MSA file is provided)

## Input Requirements

### Input Files
1. **FASTA file**: Standard protein sequence file
2. **MSA file** (optional): A3M format multiple sequence alignment

### Input Preparation
```bash
# Create FASTA file (chai-1 only recognize sequence or entity with protein, dna, rna and ligand as start)
cat > input.fasta << 'EOF'
>protein|protein_name
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLA
EOF

# Check sequence length
python -c "
from Bio import SeqIO
record = SeqIO.read('input.fasta', 'fasta')
print(f'Sequence: {record.id}')
print(f'Length: {len(record.seq)} residues')
"
```

### MSA Preparation (Optional) In general we do not generate MSA for chai-1 prediction
```bash
# Generate MSA using HHblits or similar tools
hhblits -i input.fasta -d uniprot20_2016_02 -oa3m msa.a3m -n 3

# Or use existing MSA file
cp existing_msa.a3m msa.a3m
```

## Output Structure

### Generated Files
```
output_directory/
├── pred.model_idx_0.cif      # Predicted structure model 0 (CIF format)
├── pred.model_idx_1.cif      # Predicted structure model 1
├── pred.model_idx_2.cif      # Predicted structure model 2
├── pred.model_idx_3.cif      # Predicted structure model 3
├── pred.model_idx_4.cif      # Predicted structure model 4
├── scores.model_idx_0.npz    # Scores for model 0
├── scores.model_idx_1.npz    # Scores for model 1
├── scores.model_idx_2.npz    # Scores for model 2
├── scores.model_idx_3.npz    # Scores for model 3
└── scores.model_idx_4.npz    # Scores for model 4
```

### Output Files Description

1. **CIF files** (`pred.model_idx_*.cif`):
   - Predicted protein structures in mmCIF format
   - Contains atomic coordinates and B-factors (pLDDT scores)
   - 5 models by default (0-4)

2. **NPZ files** (`scores.model_idx_*.npz`):
   - Binary numpy files containing prediction scores
   - Includes pLDDT, pTM, iPTM (if applicable) scores

## Output Interpretation

### Using Provided Utility Script
The utility script at `/data_test/home/lzzheng/apps/zDBProd/zftempl/utils.py` provides functions to analyze Chai-1 outputs:

```python
#!/usr/bin/env python
import sys
sys.path.append('/data_test/home/lzzheng/apps/zDBProd/zftempl')
from utils import calculate_confidence, load_chai_scores

# Calculate average pLDDT from CIF file
cif_file = "output_chai1/pred.model_idx_0.cif"
avg_plddt = calculate_confidence(cif_file)
print(f"Average pLDDT: {avg_plddt:.1f}")

# Load all scores from NPZ file
scores_file = "output_chai1/scores.model_idx_0.npz"
scores = load_chai_scores(scores_file)
print(f"Scores: {scores}")

# Common scores to check
if 'iptm' in scores:
    print(f"iPTM: {scores['iptm']:.3f}")
if 'ptm' in scores:
    print(f"pTM: {scores['ptm']:.3f}")
if 'aggregate_score' in scores:
    print(f"Aggregate score: {scores['aggregate_score']:.3f}")
```

### Command Line Analysis
```bash
# Analyze all models
python analyze_chai1.py --output_dir ./output_chai1 --report summary.csv

# Extract pLDDT from CIF files
for cif in output_chai1/pred.model_idx_*.cif; do
    model=$(basename "$cif" .cif | sed 's/pred.model_idx_//')
    plddt=$(python -c "
import sys
sys.path.append('/data_test/home/lzzheng/apps/zDBProd/zftempl')
from utils import calculate_confidence
print(f'{calculate_confidence(\"$cif\"):.1f}')
")
    echo "Model $model: pLDDT = $plddt"
done

# Compare model scores
python -c "
import numpy as np
import sys
sys.path.append('/data_test/home/lzzheng/apps/zDBProd/zftempl')
from utils import load_chai_scores

scores_list = []
for i in range(5):
    scores = load_chai_scores(f'output_chai1/scores.model_idx_{i}.npz')
    plddt = scores.get('plddt', 0)
    ptm = scores.get('ptm', 0)
    iptm = scores.get('iptm', 0)
    scores_list.append((i, plddt, ptm, iptm))

# Sort by pLDDT (descending)
scores_list.sort(key=lambda x: x[1], reverse=True)
print('Model ranking by pLDDT:')
for rank, (model, plddt, ptm, iptm) in enumerate(scores_list, 1):
    print(f'{rank}. Model {model}: pLDDT={plddt:.1f}, pTM={ptm:.3f}, iPTM={iptm:.3f}')
"
```

### Score Interpretation
1. **pLDDT (0-100)**:
   - >90: Very high confidence
   - 70-90: High confidence
   - 50-70: Medium confidence
   - <50: Low confidence

2. **pTM (0-1)**:
   - >0.7: Very good prediction
   - 0.5-0.7: Good prediction
   - <0.5: May have issues

3. **iPTM (0-1)**:
   - Interface predicted TM-score (for complexes)
   - >0.5: Good interface prediction

## Advanced Usage

### Batch Processing
```bash
# Process multiple FASTA files
for fasta in sequences/*.fasta; do
    base=$(basename "$fasta" .fasta)
    output_dir="predictions/${base}_chai1"
    
    # Submit to SLURM
    cd /data_test/home/lzzheng/bin
    ./submit_slurm_gpu.sh "cd /data_test/share/pub_tools/chai-lab && ./run.sh $fasta $output_dir" 4 4090
done
```

### With Custom MSA
```bash
# Generate MSA for each sequence
for fasta in sequences/*.fasta; do
    base=$(basename "$fasta" .fasta)
    
    # Generate MSA
    hhblits -i "$fasta" -d uniprot20_2016_02 -oa3m "msa/${base}.a3m" -n 3
    
    # Run Chai-1 with MSA
    cd /data_test/share/pub_tools/chai-lab
    ./run.sh "$fasta" "predictions/${base}" "msa/${base}.a3m"
done
```

### Quality Control Script
```bash
#!/bin/bash
# quality_check.sh
output_dir=$1

python -c "
import sys
sys.path.append('/data_test/home/lzzheng/apps/zDBProd/zftempl')
from utils import calculate_confidence, load_chai_scores
import numpy as np

# Check all models
for i in range(5):
    cif = f'$output_dir/pred.model_idx_{i}.cif'
    npz = f'$output_dir/scores.model_idx_{i}.npz'
    
    plddt = calculate_confidence(cif)
    scores = load_chai_scores(npz)
    
    print(f'Model {i}:')
    print(f'  pLDDT: {plddt:.1f}')
    if 'ptm' in scores:
        print(f'  pTM: {scores[\"ptm\"]:.3f}')
    if 'iptm' in scores:
        print(f'  iPTM: {scores[\"iptm\"]:.3f}')
    
    # Quality thresholds
    if plddt < 50:
        print(f'  WARNING: Low confidence (pLDDT < 50)')
    elif plddt > 70:
        print(f'  GOOD: High confidence (pLDDT > 70)')
"
```

## Troubleshooting

### Common Issues
1. **GPU memory errors**: Use smaller batch size or shorter sequences
2. **SLURM job failures**: Check partition availability and resource limits
3. **MSA server errors**: Run without MSA or provide local MSA file

### Debug Commands
```bash
# Test Chai-1 installation
cd /data_test/share/pub_tools/chai-lab
./run.sh test.fasta test_output

# Check SLURM script
cat /data_test/home/lzzheng/bin/submit_slurm_gpu.sh

# Check GPU availability
sinfo -p 4090,gpu_part

# Test utility functions
python -c "
import sys
sys.path.append('/data_test/home/lzzheng/apps/zDBProd/zftempl')
from utils import calculate_confidence
print('Utility functions loaded successfully')
"
```

## Integration with Other Tools

### Convert CIF to PDB
```bash
# Using the provided utility
python -c "
import sys
sys.path.append('/data_test/home/lzzheng/apps/zDBProd/zftempl')
from utils import convert_cif_to_pdb

cif_file = 'output_chai1/pred.model_idx_0.cif'
pdb_file = 'output_chai1/model_0.pdb'
convert_cif_to_pdb(cif_file, pdb_file)
print(f'Converted {cif_file} to {pdb_file}')
"

# Or use BioPython directly
python -c "
from Bio.PDB import MMCIFParser, PDBIO
parser = MMCIFParser()
structure = parser.get_structure('model', 'output_chai1/pred.model_idx_0.cif')
io = PDBIO()
io.set_structure(structure)
io.save('output_chai1/model_0.pdb')
"
```

### Compare with Experimental Structure
```bash
# Using TM-align
TM-align experimental.pdb output_chai1/model_0.pdb -o alignment.txt

# Or using the utility's RMSD function
python -c "
import sys
sys.path.append('/data_test/home/lzzheng/apps/zDBProd/zftempl')
from utils import calculate_rmsd_after_superimpose
from Bio.PDB import PDBParser

parser = PDBParser()
ref = parser.get_structure('ref', 'experimental.pdb')
pred = parser.get_structure('pred', 'output_chai1/model_0.pdb')

rmsd = calculate_rmsd_after_superimpose(ref, pred)
print(f'RMSD after DNA backbone superposition: {rmsd:.3f} Å')
"
```

## Performance Tips

### SLURM Optimization
```bash
# Request appropriate resources
# For small proteins (<300 residues)
./submit_slurm_gpu.sh "<command>" 2 4090

# For large proteins (>500 residues)
./submit_slurm_gpu.sh "<command>" 8 4090

# Monitor job status
squeue -u $USER
sacct -j <job_id>
```

### Memory Management
```bash
# Check memory usage
nvidia-smi

# For memory-intensive predictions
# Run with fewer parallel jobs or smaller models
```

## Python Tools

### Chai-1 Batch Processing Tool (`chai1_tool.py`)
A comprehensive Python tool for batch Chai-1 predictions is available at: `/data_test/home/lzzheng/apps/ProtDesignTools/tools/chai1_tool.py`

#### Features:
- Process single or multiple sequences from FASTA files
- Distribute jobs evenly across SLURM partitions
- Limit concurrent jobs with `--max-jobs`
- Wait for job completion and extract scores
- Generate detailed CSV reports with all 5 model scores
- Calculate average and best scores across models

#### Basic Usage:
```bash
cd /data_test/home/lzzheng/apps/ProtDesignTools
python tools/chai1_tool.py --fasta sequences.fasta --output ./results
```

#### Advanced Usage:
```bash
# Distribute across multiple partitions
python tools/chai1_tool.py --fasta sequences.fasta --output ./results \
  --partitions "4090,3090" --max-jobs 20 --ncpus 4

# With custom MSA files
python tools/chai1_tool.py --fasta sequences.fasta --output ./results \
  --msa-dir ./msa_files --msa-suffix ".a3m"

# Skip waiting for jobs
python tools/chai1_tool.py --fasta sequences.fasta --output ./results --no-wait

# Skip existing predictions
python tools/chai1_tool.py --fasta sequences.fasta --output ./results --skip-existing
```

#### Command Line Arguments:
- `--fasta`: Input FASTA file (required)
- `--output`: Output directory (default: `./chai1_results`)
- `--partitions`: Comma-separated SLURM partitions (default: "4090")
- `--max-jobs`: Maximum concurrent jobs
- `--ncpus`: CPUs per job (default: 4)
- `--msa-dir`: Directory with MSA files (optional)
- `--msa-suffix`: Suffix for MSA files (default: ".a3m")
- `--local`: Run locally instead of SLURM
- `--no-wait`: Submit jobs and exit without waiting
- `--skip-existing`: Skip sequences with existing predictions

#### Directory Structure:
```
output_directory/
├── fastas/                    # Input FASTA files (created by tool)
│   ├── sequence1.fasta
│   ├── sequence2.fasta
│   └── ...
├── sequence1/                 # Created by Chai-1
│   ├── pred.model_idx_0.cif
│   ├── pred.model_idx_1.cif
│   ├── ... (5 models total)
│   ├── scores.model_idx_0.npz
│   └── ... (5 models total)
├── sequence2/                 # Created by Chai-1
│   └── ...
└── chai1_scores_summary.csv  # Final report
```

#### Output Report:
The tool generates a CSV report (`chai1_scores_summary.csv`) with:
- For each sequence: the best model (highest pLDDT) scores only
- pLDDT, pTM, and iPTM for the best model
- Sequence information and protein sequences

## References
- Chai-1 installation: `/data_test/share/pub_tools/chai-lab/`
- SLURM submission: `/data_test/home/lzzheng/bin/submit_slurm_gpu.sh`
- Analysis utilities: `/data_test/home/lzzheng/apps/zDBProd/zftempl/utils.py`
- Chai-1 batch tool: `/data_test/home/lzzheng/apps/ProtDesignTools/tools/chai1_tool.py`
