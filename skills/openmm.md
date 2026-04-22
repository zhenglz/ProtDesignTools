# OpenMM Skill

## Overview
OpenMM is a high-performance toolkit for molecular simulation that can run on CPUs and GPUs. It's used for molecular dynamics, energy minimization, and free energy calculations.

## Installation
```bash
# Install via conda (recommended)
conda create -n openmm python=3.9
conda activate openmm
conda install -c conda-forge openmm

# Install via pip
pip install openmm

# Install with CUDA support
conda install -c conda-forge openmm cudatoolkit
```

## Basic Usage

### Command Line Interface
```bash
# Run molecular dynamics simulation
python run_simulation.py \
  --input <input_pdb> \
  --output <output_directory> \
  --length <simulation_time_ns> \
  --platform <CUDA/CPU/OpenCL>

# Example
python run_simulation.py \
  --input protein.pdb \
  --output ./simulation \
  --length 100 \
  --platform CUDA
```

### Python Script Example
```python
#!/usr/bin/env python
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit

# Load PDB
pdb = app.PDBFile('input.pdb')

# Create force field
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# Create system
system = forcefield.createSystem(pdb.topology, 
                                 nonbondedMethod=app.PME,
                                 nonbondedCutoff=1.0*unit.nanometer,
                                 constraints=app.HBonds)

# Create integrator
integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 
                                         1.0/unit.picosecond, 
                                         2.0*unit.femtoseconds)

# Create simulation
simulation = app.Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# Minimize energy
simulation.minimizeEnergy()

# Run dynamics
simulation.reporters.append(app.PDBReporter('output.pdb', 1000))
simulation.reporters.append(app.StateDataReporter('data.csv', 1000, 
                                                   step=True, time=True,
                                                   potentialEnergy=True,
                                                   temperature=True))
simulation.step(5000000)  # 10 ns
```

## Input Requirements

### Input Files
1. **PDB file**: Structure file
2. **Force field files**: XML force field definitions
3. **Parameter files**: Additional parameters if needed

### Input Preparation
```bash
# Prepare protein structure
python prepare_structure.py --pdb protein.pdb --output prepared.pdb

# Add missing residues (if needed)
python modeller_script.py --input protein.pdb --output complete.pdb

# Check protonation states
python protonate.py --input protein.pdb --output protonated.pdb
```

## Output Structure

### Generated Files
```
simulation_output/
├── trajectory.dcd                 # Trajectory file (binary)
├── trajectory.pdb                 # Trajectory file (PDB format)
├── data.csv                       # Simulation data (energy, temp, etc.)
├── restart.xml                    # Restart file
├── checkpoint.chk                 # Checkpoint file
└── log.txt                        # Simulation log
```

### Output Files Description

1. **Trajectory files**:
   - DCD: Binary format, efficient for storage
   - PDB: Text format, readable by most tools
   - Contains atomic coordinates over time

2. **Data file** (CSV):
   - Time (ps)
   - Potential energy (kJ/mol)
   - Kinetic energy (kJ/mol)
   - Temperature (K)
   - Volume (nm³)
   - Pressure (bar)

3. **Restart files**:
   - Checkpoint: Binary state for resuming
   - XML: Serialized simulation state

## Output Interpretation

### Trajectory Analysis
```bash
# Analyze RMSD
python analyze_rmsd.py --trajectory trajectory.dcd --reference reference.pdb

# Calculate radius of gyration
python analyze_rog.py --trajectory trajectory.dcd

# Extract frames
python extract_frames.py --trajectory trajectory.dcd --output frames/ --stride 100
```

### Energy Analysis
```bash
# Plot energy over time
python plot_energy.py --data data.csv --output energy_plot.png

# Check equilibration
python check_equilibration.py --data data.csv --output equilibration_report.txt
```

### Example Analysis Script
```python
import mdtraj as md
import numpy as np

# Load trajectory
traj = md.load('trajectory.dcd', top='reference.pdb')

# Calculate RMSD
rmsd = md.rmsd(traj, traj, 0)
print(f'Average RMSD: {np.mean(rmsd):.3f} nm')
print(f'RMSD fluctuation: {np.std(rmsd):.3f} nm')

# Calculate radius of gyration
rog = md.compute_rg(traj)
print(f'Average Rg: {np.mean(rog):.3f} nm')

# Calculate secondary structure
dssp = md.compute_dssp(traj)
print(f'Helix content: {(dssp == "H").mean():.2%}')
print(f'Sheet content: {(dssp == "E").mean():.2%}')
```

## Advanced Usage

### Enhanced Sampling
```bash
# Metadynamics
python run_metadynamics.py \
  --input protein.pdb \
  --cv "distance:10-20" \
  --output metadynamics

# Replica exchange
python run_replica_exchange.py \
  --input protein.pdb \
  --temperatures "300,310,320,330,340" \
  --output remd
```

### Free Energy Calculations
```bash
# Alchemical free energy
python run_free_energy.py \
  --input complex.pdb \
  --ligand ligand.pdb \
  --method "TI" \
  --output free_energy

# MM/PBSA
python run_mmpbsa.py \
  --trajectory trajectory.dcd \
  --topology topology.pdb \
  --output mmpbsa_results.csv
```

### GPU Acceleration
```bash
# Run on GPU with CUDA
python run_simulation.py --platform CUDA --device 0

# Run on multiple GPUs
python run_simulation.py --platform CUDA --devices "0,1"

# Benchmark performance
python benchmark_openmm.py --platforms "CUDA,CPU,OpenCL"
```

## Troubleshooting

### Common Issues
1. **GPU memory errors**: Reduce system size or use smaller cutoff
2. **Instability**: Reduce timestep or use constraints
3. **Force field errors**: Check atom types and parameters

### Debug Commands
```bash
# Test installation
python -c "import simtk.openmm as mm; print(mm.Platform.getPlatformByName('CUDA'))"

# Check GPU availability
python -c "
import simtk.openmm as mm
platform = mm.Platform.getPlatformByName('CUDA')
print(f'CUDA available: {platform.getNumDevices() > 0}')
"

# Validate input structure
python validate_structure.py --pdb protein.pdb
```

## Integration with Other Tools

### Analysis Pipeline
```bash
# Complete analysis workflow
python prepare_system.py --input protein.pdb --output system.pdb
python run_equilibration.py --input system.pdb --output equilibrated.pdb
python run_production.py --input equilibrated.pdb --length 100 --output trajectory.dcd
python analyze_trajectory.py --trajectory trajectory.dcd --output analysis_report.pdf
```

### Visualization
```bash
# Visualize trajectory with VMD
vmd reference.pdb trajectory.dcd

# Generate movie
python make_movie.py --trajectory trajectory.dcd --output movie.mp4

# Interactive analysis with NGLView
python interactive_analysis.py --trajectory trajectory.dcd
```

## Performance Tips

### Optimization
```bash
# Use multiple time step integrator
python run_fast.py --integrator "MTS" --timestep 4.0

# Use hydrogen mass repartitioning
python run_hmr.py --timestep 4.0

# Optimize non-bonded interactions
python run_optimized.py --cutoff 1.2 --pme
```

### Monitoring
```bash
# Monitor GPU usage
nvidia-smi -l 1

# Monitor memory usage
watch -n 1 "free -h"

# Log performance metrics
python run_with_monitoring.py --output performance.log
```

## References
- [OpenMM Documentation](http://docs.openmm.org)
- [OpenMM GitHub](https://github.com/openmm/openmm)
- [MDTraj](https://www.mdtraj.org/) for trajectory analysis