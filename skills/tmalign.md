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
TMalign <structure1> <structure2>

# Example
TMalign protein1.pdb protein2.pdb

# With output file
TMalign protein1.pdb protein2.pdb -o TM.sup

# Start with a pre-existing alignment
TMalign protein1.pdb protein2.pdb -i align.txt

# Stick to a pre-existing alignment
TMalign protein1.pdb protein2.pdb -I align.txt

# Output rotation matrix
TMalign protein1.pdb protein2.pdb -m matrix.txt

# Normalize TM-score by average length
TMalign protein1.pdb protein2.pdb -a

# Normalize TM-score by specific length (e.g., 100 residues)
TMalign protein1.pdb protein2.pdb -L 100

# Scale TM-score by specific d0 (e.g., 5 Å)
TMalign protein1.pdb protein2.pdb -d 5
```

### Key Parameters
- `-o`: Output superposition files (TM.sup, TM.sup_all, TM.sup_atm)
- `-i`: Start with alignment from file
- `-I`: Stick to alignment from file
- `-m`: Output rotation matrix
- `-a`: Normalize TM-score by average length of both proteins
- `-L`: Normalize TM-score by specified length
- `-d`: Scale TM-score by specified d0 value
- `-fast`: Fast mode (less accurate)
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
# Parse TMalign output
python -c "
import re

with open('alignment.txt', 'r') as f:
    content = f.read()
    
# Extract TM-score (two values reported)
tm_matches = re.findall(r'TM-score=\s*([\d.]+)', content)
if tm_matches:
    tm_score1 = float(tm_matches[0])  # Normalized by length of Chain_1
    tm_score2 = float(tm_matches[1])  # Normalized by length of Chain_2
    print(f'TM-score normalized by Chain_1: {tm_score1:.5f}')
    print(f'TM-score normalized by Chain_2: {tm_score2:.5f}')
    print('(Use TM-score normalized by length of the reference protein)')
    
    # Interpretation
    if tm_score1 > 0.5:
        print('Interpretation: Same fold (likely homologous)')
    elif tm_score1 > 0.3:
        print('Interpretation: Similar fold (possible homology)')
    else:
        print('Interpretation: Different fold (random structural similarity)')

# Extract RMSD
rmsd_match = re.search(r'RMSD=\s*([\d.]+)', content)
if rmsd_match:
    print(f'RMSD: {rmsd_match.group(1)} Å')

# Extract alignment length and sequence identity
len_match = re.search(r'Aligned length=\s*(\d+)', content)
if len_match:
    aligned_length = int(len_match.group(1))
    print(f'Aligned length: {aligned_length} residues')

seq_id_match = re.search(r'Seq_ID=n_identical/n_aligned=\s*([\d.]+)', content)
if seq_id_match:
    seq_id = float(seq_id_match.group(1))
    print(f'Sequence identity: {seq_id:.3f}')

# Extract chain lengths
len1_match = re.search(r'Length of Chain_1:\s*(\d+)', content)
len2_match = re.search(r'Length of Chain_2:\s*(\d+)', content)
if len1_match and len2_match:
    len1 = int(len1_match.group(1))
    len2 = int(len2_match.group(1))
    print(f'Chain_1 length: {len1} residues')
    print(f'Chain_2 length: {len2} residues')
    print(f'Alignment coverage (Chain_1): {aligned_length/len1*100:.1f}%')
    print(f'Alignment coverage (Chain_2): {aligned_length/len2*100:.1f}%')
"
```

### Visualizing Alignment
```bash
# View superimposed structures in PyMOL
pymol TM.sup_all

# Or using the generated PyMOL scripts
pymol -d @TM.sup.pml              # View aligned regions (C-alpha)
pymol -d @TM.sup_all.pml          # View all regions (C-alpha)
pymol -d @TM.sup_atm.pml          # View aligned regions (full-atom)
pymol -d @TM.sup_all_atm.pml      # View all regions (full-atom)
pymol -d @TM.sup_all_atm_lig.pml  # View all regions with ligands

# Or generate alignment figure
python visualize_alignment.py \
  --structure1 TM.sup \
  --structure2 protein2.pdb \
  --output alignment.png
```

## Per-Residue Pairing Analysis

### Understanding the Alignment Output
TMalign provides a visual alignment string showing residue pairs:
- `:` denotes aligned residue pairs with distance < 5.0 Å
- `.` denotes other aligned residues
- Each character corresponds to one residue in the alignment

### Extracting Per-Residue Correspondence
```python
import re
from collections import defaultdict

def parse_tmalign_output(tmalign_output):
    """Parse TMalign output to extract per-residue alignment information."""
    
    lines = tmalign_output.strip().split('\n')
    
    # Find the alignment section
    alignment_start = None
    for i, line in enumerate(lines):
        if line.startswith('(') and ':' in line and '.' in line:
            alignment_start = i
            break
    
    if alignment_start is None:
        raise ValueError("Could not find alignment section in TMalign output")
    
    # Extract sequences and alignment string
    seq1_line = lines[alignment_start - 2]
    seq2_line = lines[alignment_start - 1]
    align_line = lines[alignment_start]
    
    # Remove the prefix before the actual alignment
    seq1 = seq1_line.split()[-1] if ' ' in seq1_line else seq1_line
    seq2 = seq2_line.split()[-1] if ' ' in seq2_line else seq2_line
    alignment = align_line.split()[-1] if ' ' in align_line else align_line
    
    # Parse per-residue correspondence
    residue_pairs = []
    seq1_idx = 0
    seq2_idx = 0
    
    for align_char in alignment:
        if align_char == ':' or align_char == '.':
            # This is an aligned position
            if seq1_idx < len(seq1) and seq2_idx < len(seq2):
                residue_pairs.append({
                    'chain1_residue': seq1[seq1_idx],
                    'chain1_position': seq1_idx + 1,
                    'chain2_residue': seq2[seq2_idx],
                    'chain2_position': seq2_idx + 1,
                    'alignment_type': 'close' if align_char == ':' else 'distant',
                    'distance_threshold': '<5.0A' if align_char == ':' else '>=5.0A'
                })
            seq1_idx += 1
            seq2_idx += 1
        elif align_char == '-':
            # Gap in sequence 1
            seq2_idx += 1
        else:
            # Gap in sequence 2 or other character
            seq1_idx += 1
    
    return {
        'sequence1': seq1,
        'sequence2': seq2,
        'alignment_string': alignment,
        'residue_pairs': residue_pairs,
        'num_aligned': len([c for c in alignment if c in ':.']),
        'num_close_pairs': alignment.count(':'),
        'num_distant_pairs': alignment.count('.')
    }

def get_residue_mapping(tmalign_output, reference_chain=1):
    """Get residue index mapping between two structures.
    
    Args:
        tmalign_output: Raw TMalign output text
        reference_chain: Which chain to use as reference (1 or 2)
    
    Returns:
        dict: Mapping from reference residue index to target residue index
    """
    alignment_info = parse_tmalign_output(tmalign_output)
    
    mapping = {}
    for pair in alignment_info['residue_pairs']:
        if reference_chain == 1:
            mapping[pair['chain1_position']] = pair['chain2_position']
        else:
            mapping[pair['chain2_position']] = pair['chain1_position']
    
    return mapping

def calculate_per_residue_rmsd(pdb1, pdb2, alignment_info):
    """Calculate per-residue RMSD based on TMalign alignment."""
    from Bio.PDB import PDBParser, Superimposer
    import numpy as np
    
    parser = PDBParser()
    struct1 = parser.get_structure('chain1', pdb1)
    struct2 = parser.get_structure('chain2', pdb2)
    
    # Get C-alpha atoms
    ca_atoms1 = [atom for atom in struct1.get_atoms() if atom.get_name() == 'CA']
    ca_atoms2 = [atom for atom in struct2.get_atoms() if atom.get_name() == 'CA']
    
    per_residue_rmsd = []
    
    for pair in alignment_info['residue_pairs']:
        idx1 = pair['chain1_position'] - 1  # Convert to 0-based
        idx2 = pair['chain2_position'] - 1
        
        if idx1 < len(ca_atoms1) and idx2 < len(ca_atoms2):
            coord1 = ca_atoms1[idx1].get_coord()
            coord2 = ca_atoms2[idx2].get_coord()
            distance = np.linalg.norm(coord1 - coord2)
            per_residue_rmsd.append({
                'chain1_residue': pair['chain1_residue'],
                'chain1_position': pair['chain1_position'],
                'chain2_residue': pair['chain2_residue'],
                'chain2_position': pair['chain2_position'],
                'distance': distance,
                'alignment_type': pair['alignment_type']
            })
    
    return per_residue_rmsd

# Example usage
if __name__ == "__main__":
    # Run TMalign and capture output
    import subprocess
    result = subprocess.run(
        ['TMalign', 'chai1_RpNMT.pdb', 'receptor.pdb'],
        capture_output=True,
        text=True
    )
    
    # Parse the output
    alignment_info = parse_tmalign_output(result.stdout)
    
    print(f"Total aligned residues: {alignment_info['num_aligned']}")
    print(f"Close pairs (<5Å): {alignment_info['num_close_pairs']}")
    print(f"Distant pairs (≥5Å): {alignment_info['num_distant_pairs']}")
    
    # Get residue mapping
    mapping = get_residue_mapping(result.stdout, reference_chain=1)
    print(f"\nResidue mapping (Chain1 → Chain2):")
    for ref_pos, target_pos in list(mapping.items())[:10]:  # Show first 10
        print(f"  {ref_pos} → {target_pos}")
    
    # Calculate per-residue distances
    per_residue_dist = calculate_per_residue_rmsd(
        'chai1_RpNMT.pdb', 
        'receptor.pdb', 
        alignment_info
    )
    
    print(f"\nPer-residue distances (first 10):")
    for dist_info in per_residue_dist[:10]:
        print(f"  {dist_info['chain1_residue']}{dist_info['chain1_position']} ↔ "
              f"{dist_info['chain2_residue']}{dist_info['chain2_position']}: "
              f"{dist_info['distance']:.2f}Å ({dist_info['alignment_type']})")
```

### Using the Rotation Matrix for Structural Superposition
```python
def apply_rotation_matrix(pdb_file, matrix_file, output_file):
    """Apply TMalign rotation matrix to a structure."""
    import numpy as np
    from Bio.PDB import PDBParser, PDBIO
    
    # Parse rotation matrix
    with open(matrix_file, 'r') as f:
        lines = f.readlines()
    
    # TMalign matrix format: 3x3 rotation matrix + translation vector
    rotation = np.zeros((3, 3))
    translation = np.zeros(3)
    
    for i in range(3):
        parts = lines[i].strip().split()
        rotation[i] = [float(x) for x in parts[:3]]
        translation[i] = float(parts[3])
    
    # Parse and transform structure
    parser = PDBParser()
    structure = parser.get_structure('input', pdb_file)
    
    for atom in structure.get_atoms():
        coord = np.array(atom.get_coord())
        new_coord = np.dot(rotation, coord) + translation
        atom.set_coord(new_coord)
    
    # Save transformed structure
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)
    
    return output_file

# Example usage
apply_rotation_matrix('chai1_RpNMT.pdb', 'matrix.txt', 'chai1_aligned.pdb')
```

## Practical Examples and Interpretation

### Example 1: High Similarity (Same Protein)
```bash
TMalign chai1_RpNMT.pdb receptor.pdb
```
**Output Interpretation:**
- TM-score: 0.96491 (very high, >0.5 indicates same fold)
- RMSD: 1.22 Å (low, indicates good structural match)
- Aligned length: 286 residues (100% coverage)
- Sequence identity: 1.000 (identical sequences)
- **Conclusion**: These are essentially the same structure with minor variations

### Example 2: Moderate Similarity (Related Proteins)
```bash
TMalign proteinA.pdb proteinB.pdb
```
**Typical Output:**
- TM-score: 0.45 (moderate, 0.3-0.5 range)
- RMSD: 3.5 Å (moderate)
- Aligned length: 150 residues (60% coverage)
- Sequence identity: 0.25 (low sequence similarity)
- **Conclusion**: Possibly homologous proteins with similar fold but divergent sequences

### Example 3: Low Similarity (Different Folds)
```bash
TMalign alpha_helix_bundle.pdb beta_barrel.pdb
```
**Typical Output:**
- TM-score: 0.20 (low, <0.3)
- RMSD: 8.2 Å (high)
- Aligned length: 80 residues (30% coverage)
- Sequence identity: 0.10 (very low)
- **Conclusion**: Different structural folds

## Advanced Usage

### Multiple Structure Alignment
```bash
# Align multiple structures to a reference
for pdb in structures/*.pdb; do
    base=$(basename "$pdb" .pdb)
    TMalign reference.pdb "$pdb" -o "alignments/${base}_alignment.txt"
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

## Practical Examples and Interpretation

### Example 1: High Similarity (Same Protein)
```bash
TMalign chai1_RpNMT.pdb receptor.pdb
```
**Output Interpretation:**
- TM-score: 0.96491 (very high, >0.5 indicates same fold)
- RMSD: 1.22 Å (low, indicates good structural match)
- Aligned length: 286 residues (100% coverage)
- Sequence identity: 1.000 (identical sequences)
- **Conclusion**: These are essentially the same structure with minor variations

### Example 2: Moderate Similarity (Related Proteins)
```bash
TMalign proteinA.pdb proteinB.pdb
```
**Typical Output:**
- TM-score: 0.45 (moderate, 0.3-0.5 range)
- RMSD: 3.5 Å (moderate)
- Aligned length: 150 residues (60% coverage)
- Sequence identity: 0.25 (low sequence similarity)
- **Conclusion**: Possibly homologous proteins with similar fold but divergent sequences

### Example 3: Low Similarity (Different Folds)
```bash
TMalign alpha_helix_bundle.pdb beta_barrel.pdb
```
**Typical Output:**
- TM-score: 0.20 (low, <0.3)
- RMSD: 8.2 Å (high)
- Aligned length: 80 residues (30% coverage)
- Sequence identity: 0.10 (very low)
- **Conclusion**: Different structural folds

## Common Use Cases

### 1. Structure Comparison for Homology Detection
```bash
# Compare novel structure to known structures
TMalign novel_protein.pdb known_template.pdb -o comparison.txt

# Check if TM-score > 0.5 (likely homologous)
python -c "
import re
with open('comparison.txt') as f:
    content = f.read()
tm_match = re.search(r'TM-score=\s*([\d.]+)', content)
if tm_match and float(tm_match.group(1)) > 0.5:
    print('Likely homologous structures')
else:
    print('Different folds')
"
```

### 2. Structural Superposition for Visualization
```bash
# Generate superposition for PyMOL
TMalign mobile.pdb target.pdb -o superposed.pdb

# View in PyMOL
pymol superposed.pdb target.pdb
```

### 3. Residue Correspondence for Mutational Analysis
```python
# Map residues between wild-type and mutant structures
alignment = parse_tmalign_output(tmalign_output)
for pair in alignment['residue_pairs']:
    if pair['chain1_residue'] != pair['chain2_residue']:
        print(f"Mutation: {pair['chain1_residue']}{pair['chain1_position']}"
              f" → {pair['chain2_residue']}{pair['chain2_position']}")
```

### 4. Quality Assessment of Structure Predictions
```bash
# Compare predicted model to experimental structure
TMalign predicted.pdb experimental.pdb

# Good prediction: TM-score > 0.8, RMSD < 2.0 Å
# Acceptable prediction: TM-score > 0.5, RMSD < 4.0 Å
```

## Quick Reference

### TM-score Interpretation
- **> 0.5**: Same fold (likely homologous)
- **0.3-0.5**: Similar fold (possible homology)
- **< 0.3**: Different fold (random structural similarity)

### Key Output Files
- `TM.sup`: Superposed C-alpha traces (aligned regions)
- `TM.sup_all`: Superposed C-alpha traces (all regions)
- `TM.sup_atm`: Superposed full-atom structures (aligned regions)
- `matrix.txt`: Rotation/translation matrix (with `-m` flag)

### Useful One-Liners
```bash
# Get just the TM-score
TMalign struct1.pdb struct2.pdb 2>&1 | grep "TM-score" | head -1

# Batch process multiple comparisons
for i in *.pdb; do for j in *.pdb; do [ "$i" != "$j" ] && echo -n "$i vs $j: " && TMalign "$i" "$j" 2>&1 | grep "TM-score" | head -1; done; done

# Create alignment matrix for multiple structures
python -c "
import itertools, subprocess, re
pdbs = ['struct1.pdb', 'struct2.pdb', 'struct3.pdb']
for a,b in itertools.combinations(pdbs, 2):
    out = subprocess.run(['TMalign', a, b], capture_output=True, text=True)
    tm = re.search(r'TM-score=\\s*([\\d.]+)', out.stdout)
    print(f'{a} vs {b}: {tm.group(1) if tm else \"N/A\"}')
"
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