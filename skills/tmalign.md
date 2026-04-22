# TM-align Skill

## Overview
TM-align is a protein structure alignment algorithm that uses TM-score to evaluate structural similarity. It's widely used for comparing protein structures and detecting structural homology.

## Installation
```bash
# Download TM-align binary
wget https://zhanggroup.org/TM-align/TM-align.cpp
g++ -static -O3 -ffast-math -lm -o TM-align TM-align.cpp

# Or download pre-compiled binary
wget https://zhanggroup.org/TM-align/TM-align.gz
gunzip TM-align.gz
chmod +x TM-align
sudo mv TM-align /usr/local/bin/

# Python wrapper (optional)
pip install tmscoring
```

## Basic Usage

### Command Line Interface
```bash
# Align two structures
TM-align <structure1> <structure2>

# Example
TM-align protein1.pdb protein2.pdb

# With output file
TM-align protein1.pdb protein2.pdb -o alignment.txt
```

### Key Parameters
- `-o`: Output file for alignment details
- `-a`: Output aligned structures (TM.sup, TM.sup_all)
- `-d`: Output distance matrix
- `-m`: Output rotation matrix
- `-fast`: Fast mode (less accurate)
- `-L`: Length of structure to align (default: length of shorter structure)
- `-seqfix`: Fix sequence order (for comparing identical sequences)

## Input Requirements

### Input Files
1. **PDB files**: Standard protein structure files (.pdb format)
2. **Can also accept**: mmCIF, PDBx files

### Input Preparation
```bash
# Clean PDB files (remove heteroatoms, waters)
grep -E "^ATOM" protein1.pdb > protein1_clean.pdb
grep -E "^ATOM" protein2.pdb > protein2_clean.pdb

# Extract specific chains (if needed)
grep -E "^ATOM.*A" protein1.pdb > protein1_chainA.pdb

# Check structure integrity
python -c "
from Bio.PDB import PDBParser
parser = PDBParser()
s1 = parser.get_structure('s1', 'protein1.pdb')
s2 = parser.get_structure('s2', 'protein2.pdb')
print(f'Structure 1: {len(list(s1.get_residues()))} residues')
print(f'Structure 2: {len(list(s2.get_residues()))} residues')
"
```

## Output Structure

### Generated Files
```
alignment_results/
├── TM.sup                    # Superimposed structure 1
├── TM.sup_all               # Superimposed structure 2
├── alignment.txt            # Detailed alignment output
├── rotation_matrix.txt      # Rotation matrix (if -m specified)
└── distance_matrix.txt      # Distance matrix (if -d specified)
```

### Output Files Description

1. **TM.sup files**:
   - Superimposed structures in PDB format
   - TM.sup: Structure 1 rotated to match structure 2
   - TM.sup_all: Both structures superimposed

2. **Alignment file**:
   - TM-score and RMSD values
   - Alignment length and coverage
   - Sequence alignment
   - Structural alignment details

## Output Interpretation

### Key Metrics
1. **TM-score (0-1)**:
   - >0.5: Same fold (likely homologous)
   - 0.3-0.5: Similar fold (possible homology)
   - <0.3: Different fold

2. **RMSD (Å)**:
   - Lower values indicate better alignment
   - Normalized by alignment length

3. **Alignment coverage**:
   - Percentage of residues aligned
   - IDEN: Sequence identity of aligned residues

### Example Output Analysis
```bash
# Parse TM-align output
python -c "
import re

with open('alignment.txt', 'r') as f:
    content = f.read()
    
# Extract TM-score
tm_match = re.search(r'TM-score\s*=\s*([\d.]+)', content)
if tm_match:
    tm_score = float(tm_match.group(1))
    print(f'TM-score: {tm_score:.3f}')
    if tm_score > 0.5:
        print('Interpretation: Same fold (likely homologous)')
    elif tm_score > 0.3:
        print('Interpretation: Similar fold (possible homology)')
    else:
        print('Interpretation: Different fold')

# Extract RMSD
rmsd_match = re.search(r'RMSD\s*=\s*([\d.]+)', content)
if rmsd_match:
    print(f'RMSD: {rmsd_match.group(1)} Å')

# Extract alignment length
len_match = re.search(r'Aligned length=\s*(\d+)', content)
if len_match:
    print(f'Aligned length: {len_match.group(1)} residues')
"
```

### Visualizing Alignment
```bash
# View superimposed structures in PyMOL
pymol TM.sup_all

# Or generate alignment figure
python visualize_alignment.py \
  --structure1 TM.sup \
  --structure2 protein2.pdb \
  --output alignment.png
```

## Advanced Usage

### Multiple Structure Alignment
```bash
# Align multiple structures to a reference
for pdb in structures/*.pdb; do
    base=$(basename "$pdb" .pdb)
    TM-align reference.pdb "$pdb" -o "alignments/${base}_alignment.txt"
done

# Analyze all alignments
python analyze_multiple_alignments.py alignments/
```

### Chain-specific Alignment
```bash
# Align specific chains
TM-align protein1_chainA.pdb protein2_chainB.pdb

# Or extract chains on the fly
TM-align <(grep -E "^ATOM.*A" protein1.pdb) <(grep -E "^ATOM.*B" protein2.pdb)
```

### Custom Length Normalization
```bash
# Use specific length for TM-score normalization
TM-align protein1.pdb protein2.pdb -L 100

# Compare different normalization methods
python compare_tm_scores.py --structures "protein1.pdb,protein2.pdb,protein3.pdb"
```

## Troubleshooting

### Common Issues
1. **No atoms in input**: Check PDB file format
2. **Memory errors**: Use smaller structures or `-fast` mode
3. **Incorrect chain selection**: Clean PDB files before alignment

### Debug Commands
```bash
# Test with simple structures
TM-align test1.pdb test2.pdb

# Check PDB file contents
head -20 protein1.pdb

# Validate structure
python validate_pdb.py --pdb protein1.pdb
```

## Integration with Other Tools

### Structural Classification
```bash
# Classify structures by fold
python classify_folds.py \
  --query query.pdb \
  --database fold_database/ \
  --output classification.csv

# Build dendrogram of structural similarity
python build_dendrogram.py --alignments alignments/ --output dendrogram.png
```

### Quality Assessment
```bash
# Compare model to native structure
TM-align model.pdb native.pdb -o model_vs_native.txt

# Calculate GDT scores
python calculate_gdt.py --alignment model_vs_native.txt

# Generate quality report
python generate_quality_report.py \
  --model model.pdb \
  --native native.pdb \
  --output quality_report.pdf
```

### Database Searching
```bash
# Search structural database
python search_structural_db.py \
  --query query.pdb \
  --database pdb70 \
  --output hits.csv

# Filter by TM-score
python filter_by_tmscore.py --hits hits.csv --threshold 0.5 --output filtered_hits.csv
```

## Performance Tips

### Batch Processing
```bash
# Parallel alignment
parallel -j 4 TM-align reference.pdb {} -o alignments/{/.}.txt ::: structures/*.pdb

# Use fast mode for screening
TM-align protein1.pdb protein2.pdb -fast
```

### Memory Management
```bash
# For very large structures
TM-align large1.pdb large2.pdb -fast -L 500

# Split large structures
python split_structure.py --pdb large.pdb --output chunks/ --size 500
```

## Alternative Tools

### Related Software
```bash
# FATCAT (flexible structure alignment)
FATCAT -p1 protein1.pdb -p2 protein2.pdb

# CE (Combinatorial Extension)
ce protein1.pdb protein2.pdb

# DALI (distance matrix alignment)
dali protein1.pdb protein2.pdb
```

## References
- [TM-align Website](https://zhanggroup.org/TM-align/)
- [TM-score Paper](https://doi.org/10.1093/nar/gki524)
- [TM-align GitHub](https://github.com/kad-ecoli/TM-align)