# RFDiffusion Skill

## Overview
RFDiffusion is a diffusion model for protein structure generation that can create novel protein structures from scratch or conditionally modify existing structures.

## Installation
```bash
# Clone the repository
git clone https://github.com/RosettaCommons/RFdiffusion.git
cd RFdiffusion

# Install dependencies (requires Conda)
conda env create -f environment.yml
conda activate rfdiffusion

# Download model weights
python scripts/download_models.py
```

## Basic Usage

### Command Line Interface
```bash
# Generate novel protein structures
python scripts/run_inference.py \
  --contigs <contig_description> \
  --out <output_prefix> \
  --num_designs <number_of_designs>

# Example: generate a 100-residue protein
python scripts/run_inference.py \
  --contigs "100" \
  --out ./output/novel_protein \
  --num_designs 10
```

### Key Parameters
- `--contigs`: Contig description string (e.g., "80-100" or "A1-100,B101-150")
- `--out`: Output file prefix
- `--num_designs`: Number of structures to generate
- `--inpaint_str`: Residues to redesign (e.g., "A10-30")
- `--hotspots`: Residues to bias sampling around
- `--symmetry`: Symmetry specification (e.g., "C3" for cyclic symmetry)
- `--guide_scale`: Guidance scale for conditional generation

## Input Requirements

### Contig Notation
- `"100"`: Generate a 100-residue protein
- `"80-100"`: Generate protein between 80-100 residues
- `"A1-100"`: Generate chain A with residues 1-100
- `"A1-100,B101-150"`: Generate two chains
- `"A10-30/A40-60"`: Generate two segments with gap between

### Input Preparation
```bash
# For conditioning on existing structure
python scripts/run_inference.py \
  --contigs "A1-100" \
  --inpaint_str "A30-50" \
  --pdb <input_pdb> \
  --out ./output/redesigned
```

## Output Structure

### Generated Files
```
output_directory/
├── <prefix>_<design_number>.pdb          # Generated PDB structures
├── <prefix>_<design_number>.npz          # Raw model outputs
├── <prefix>_scores.csv                   # Design scores
└── <prefix>_log.txt                      # Run log
```

### Output Files Description

1. **PDB files** (`*.pdb`):
   - Generated protein structures
   - Contains atomic coordinates

2. **NPZ files** (`*.npz`):
   - Raw model outputs (logits, features)
   - Useful for debugging and analysis

3. **CSV scores** (`*_scores.csv`):
   - Contains per-design scores:
     - `plddt`: Per-residue confidence score (0-100)
     - `ptm`: Predicted TM-score
     - `iptm`: Interface predicted TM-score (for complexes)
     - `rmsd`: RMSD to input (if applicable)

## Output Interpretation

### Quality Metrics
1. **pLDDT (0-100)**:
   - >90: Very high confidence
   - 70-90: High confidence
   - 50-70: Medium confidence
   - <50: Low confidence

2. **Predicted TM-score (0-1)**:
   - >0.5: Likely correct fold
   - <0.5: May be incorrect fold

3. **Interface pTM (for complexes)**:
   - >0.5: Good interface prediction

### Example Output Analysis
```bash
# Analyze pLDDT scores
python -c "
import numpy as np
scores = np.load('output/design_0.npz')
plddt = scores['plddt']
print(f'Mean pLDDT: {plddt.mean():.1f}')
print(f'Min pLDDT: {plddt.min():.1f}')
print(f'Max pLDDT: {plddt.max():.1f}')
print(f'Residues with pLDDT < 50: {(plddt < 50).sum()}')
"

# Visualize confidence
python scripts/visualize_plddt.py --pdb output/design_0.pdb --out confidence.png
```

## Advanced Usage

### Symmetric Oligomers
```bash
# Generate C3 symmetric trimer
python scripts/run_inference.py \
  --contigs "100" \
  --symmetry "C3" \
  --out ./output/trimer \
  --num_designs 5
```

### Binding Site Design
```bash
# Design binding site around specified residues
python scripts/run_inference.py \
  --contigs "A1-100" \
  --hotspots "A30,A35,A40" \
  --out ./output/binding_site \
  --num_designs 20
```

### Scaffold from Motif
```bash
# Scaffold a functional motif
python scripts/run_inference.py \
  --contigs "A1-50/10/A51-100" \
  --inpaint_str "A1-50/A51-100" \
  --pdb motif.pdb \
  --out ./output/scaffolded_motif
```

## Troubleshooting

### Common Issues
1. **Memory errors**: Reduce `--num_designs` or use smaller contigs
2. **Contig parsing errors**: Check contig string format
3. **Model loading errors**: Ensure model weights are downloaded

### Debug Commands
```bash
# Test with small design
python scripts/run_inference.py --contigs "50" --num_designs 1 --out test

# Check PDB file
python -c "
from Bio.PDB import PDBParser
parser = PDBParser()
structure = parser.get_structure('test', 'input.pdb')
print(f'Number of chains: {len(list(structure.get_chains()))}')
print(f'Number of residues: {len(list(structure.get_residues()))}')
"
```

## Integration with Other Tools

### Downstream Analysis
```bash
# Refine with Rosetta
for pdb in output/*.pdb; do
    rosetta_scripts @flags -s "$pdb" -out:prefix "${pdb%.pdb}_refined"
done

# Score with AlphaFold2
python run_alphafold.py --pdb_dir ./output --output_dir ./alphafold_scores
```

## References
- [RFDiffusion GitHub](https://github.com/RosettaCommons/RFdiffusion)
- [RFDiffusion Paper](https://www.science.org/doi/10.1126/science.adf6591)