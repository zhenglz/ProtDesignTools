# AutoDock Vina Skill

## Overview
AutoDock Vina is a molecular docking program for predicting how small molecules (ligands) bind to protein targets. It's widely used for virtual screening and drug discovery.

## Installation
```bash
# Download and install AutoDock Vina
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
chmod +x vina_1.2.5_linux_x86_64
sudo mv vina_1.2.5_linux_x86_64 /usr/local/bin/vina

# Install Python wrapper (optional)
pip install vina
```

## Basic Usage

### Command Line Interface
```bash
# Basic docking
vina \
  --receptor <receptor_pdbqt> \
  --ligand <ligand_pdbqt> \
  --config <config_file> \
  --out <output_pdbqt> \
  --log <log_file>

# Example
vina \
  --receptor protein.pdbqt \
  --ligand ligand.pdbqt \
  --config config.txt \
  --out docked_ligand.pdbqt \
  --log docking.log
```

### Key Parameters
- `--receptor`: Receptor PDBQT file
- `--ligand`: Ligand PDBQT file
- `--config`: Configuration file with docking box parameters
- `--out`: Output PDBQT file
- `--log`: Log file
- `--cpu`: Number of CPUs to use
- `--exhaustiveness`: Search exhaustiveness (higher = more thorough)
- `--num_modes`: Number of binding modes to output

## Input Requirements

### Input File Preparation

#### 1. Prepare Receptor
```bash
# Convert PDB to PDBQT
prepare_receptor4.py -r protein.pdb -o protein.pdbqt

# Add hydrogens and charges
pythonsh prepare_receptor4.py -r protein.pdb -o protein.pdbqt -A checkhydrogens
```

#### 2. Prepare Ligand
```bash
# Convert ligand to PDBQT
obabel ligand.sdf -O ligand.pdbqt -h --partialcharge gasteiger

# Or use AutoDockTools
prepare_ligand4.py -l ligand.pdb -o ligand.pdbqt
```

#### 3. Create Configuration File
```bash
# Generate config file with docking box
cat > config.txt << 'EOF'
center_x = 10.0
center_y = 20.0
center_z = 30.0
size_x = 25.0
size_y = 25.0
size_z = 25.0
exhaustiveness = 8
num_modes = 9
energy_range = 4
EOF

# Or use AutoDockTools to define box
write_gpf4.py -l ligand.pdbqt -r protein.pdbqt
```

## Output Structure

### Generated Files
```
docking_results/
├── docked_ligand.pdbqt          # Docked ligand poses
├── docking.log                  # Detailed docking log
└── (optional) separated poses as individual files
```

### Output Files Description

1. **PDBQT output file**:
   - Contains multiple docking poses (MODEL/ENDMDL records)
   - Each pose includes estimated free energy (kcal/mol)
   - RMSD values relative to best pose

2. **Log file**:
   - Docking scores for each pose
   - RMSD table between poses
   - Runtime information

## Output Interpretation

### Analyzing Results
```bash
# Extract docking scores
grep -A 10 "mode |   affinity" docking.log

# Parse output with Python
python -c "
import re

with open('docking.log', 'r') as f:
    content = f.read()
    
# Extract scores
pattern = r'mode\s+\|\s+affinity\s+\|\s+rmsd l.b.\s+\|\s+rmsd u.b.'
matches = re.findall(r'(\d+)\s+\|\s+([-\d.]+)\s+\|\s+([\d.]+)\s+\|\s+([\d.]+)', content)
for mode, affinity, rmsd_lb, rmsd_ub in matches:
    print(f'Mode {mode}: Affinity = {affinity} kcal/mol, RMSD = {rmsd_lb}-{rmsd_ub} Å')
"

# Convert PDBQT to PDB for visualization
obabel docked_ligand.pdbqt -O docked_poses.pdb -m
```

### Score Interpretation
- **Affinity (kcal/mol)**: More negative = stronger binding
- **Typical ranges**:
  - -12 to -15: Very strong binder
  - -9 to -12: Strong binder  
  - -6 to -9: Moderate binder
  - -3 to -6: Weak binder
  - > -3: Very weak or non-binder

## Advanced Usage

### Virtual Screening
```bash
# Screen multiple ligands
for ligand in ligands/*.pdbqt; do
    base=$(basename "$ligand" .pdbqt)
    vina \
      --receptor protein.pdbqt \
      --ligand "$ligand" \
      --config config.txt \
      --out "results/${base}_docked.pdbqt" \
      --log "results/${base}.log"
done

# Analyze screening results
python analyze_screening.py results/
```

### Flexible Residue Docking
```bash
# Prepare flexible receptor
prepare_flexreceptor4.py -r protein.pdbqt -s flexible_residues.txt -g protein_flex.gpf

# Run docking with flexibility
vina \
  --receptor protein_rigid.pdbqt \
  --flex protein_flex.pdbqt \
  --ligand ligand.pdbqt \
  --config config.txt
```

### Custom Scoring Function
```bash
# Use custom scoring function weights
cat > custom_weights.txt << 'EOF'
gauss1     -0.035579
gauss2     -0.005156
repulsion   0.840245
hydrophobic -0.035069
hydrogen    -0.587439
rot         0.05846
EOF

vina --scoring custom_weights.txt --receptor protein.pdbqt --ligand ligand.pdbqt
```

## Troubleshooting

### Common Issues
1. **PDBQT conversion errors**: Check input file format and hydrogens
2. **Docking box too small**: Increase size_x, size_y, size_z
3. **No poses found**: Increase exhaustiveness or adjust box center

### Debug Commands
```bash
# Check PDBQT file
head -20 ligand.pdbqt

# Test with simple ligand
vina --receptor protein.pdbqt --ligand test_ligand.pdbqt --config test_config.txt

# Verify box parameters
python -c "
# Visualize docking box
import nglview as nv
import pytraj as pt

traj = pt.load('protein.pdb')
view = nv.show_pytraj(traj)
view.add_representation('licorice', selection='protein')
view.center()
view
"
```

## Integration with Other Tools

### Post-Docking Analysis
```bash
# Calculate binding energy with MM/PBSA
python mm_pbsa.py --complex complex.pdb --receptor protein.pdb --ligand ligand.pdb

# Visualize poses
pymol protein.pdb docked_poses.pdb

# Generate interaction diagrams
python scripts/generate_interactions.py --pose docked_pose1.pdb
```

### Workflow Automation
```bash
#!/bin/bash
# Automated docking pipeline
prepare_receptor4.py -r "$1" -o receptor.pdbqt
prepare_ligand4.py -l "$2" -o ligand.pdbqt

# Auto-detect binding site
python detect_binding_site.py --receptor receptor.pdbqt --ligand ligand.pdbqt > config.txt

vina --receptor receptor.pdbqt --ligand ligand.pdbqt --config config.txt --out docked.pdbqt

# Analyze results
python analyze_docking.py --log docking.log --output summary.csv
```

## References
- [AutoDock Vina GitHub](https://github.com/ccsb-scripps/AutoDock-Vina)
- [AutoDockTools](https://autodocksuite.scripps.edu/)
- [Open Babel](http://openbabel.org/)