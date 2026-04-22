# FoldX Scorer Skill

## Overview
FoldX is a widely used computational tool for predicting the stability changes of proteins upon mutation (ΔΔG). This skill provides comprehensive guidance on using the FoldX implementation in the zMutScan framework for mutation stability prediction.

## Key Features

- **ΔΔG prediction**: Calculates free energy change upon mutation
- **Force field-based**: Empirical force field optimized on experimental data
- **Structure-based**: Requires 3D protein structure
- **Multi-chain support**: Handles homodimers and complexes
- **Experimental validation**: Well-correlated with experimental ΔΔG values

## Prerequisites

### FoldX Installation
```bash
# Check if FoldX is installed
which FoldX  # or which foldx

# Expected location: /path/to/foldx/bin/FoldX
# If not installed, download from: http://foldxsuite.org/
```

### Required Files
1. **FoldX binary**: `FoldX` or `foldx`
2. **Rotabase library**: `rotabase.txt` (usually in same directory as FoldX)
3. **PDB file**: Clean protein structure file

## Basic Usage

### Command Line Interface
```bash
# Basic usage
/data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/foldx_scorer.py \
    <pdb_file> \
    <mutation> \
    <output_dir>

# Example
/data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/foldx_scorer.py \
    protein.pdb \
    M198F \
    ./foldx_results

# With multiple runs (default: 5)
/data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/foldx_scorer.py \
    --num_runs 10 \
    protein.pdb \
    M198F \
    ./foldx_results

# For homodimers
/data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/foldx_scorer.py \
    --is_homodimer \
    protein.pdb \
    M198F \
    ./foldx_results
```

### Python API
```python
#!/usr/bin/env python3
from zMutScan.scripts.foldx_scorer import FoldXScorer
from zMutScan.scripts.base_scorer import Mutation

# Initialize scorer
scorer = FoldXScorer(
    config={
        "args": {
            "num_runs": 5,  # Number of independent runs
            "rotabase_location": "/path/to/rotabase.txt"  # Optional
        }
    },
    output_dir="./foldx_results",
    wt_sequence="MRSLTAKEYVEAFKSFLDHSTEHQCMEAFNKQEMPHIMAGLGNGKSTLNVLGVGSGTGEQDLKMIQILQAEHPGVFINAEIIEPNPQHVAAYKELVNQAPGLQNVSFIWHQLTSSEYEQQMKEKSTHKKFDFIHMIQMLYRVEDIPNTIKFFHSCLDHHGKLLIIILSESSGWATLWKKYRHCLPLTDSGHYITSNDIEDILKRIGVEYHVFEFPSGWDITECFIEGDLVGGHMMDFLTGTKNFLGTAPVDLKKRLQEALCQPECSSRKDGRIIFNNNLSMIVVES",
    pdb_file="protein.pdb",
    is_homodimer=False  # Set to True for homodimers
)

# Score single mutation
mutation = Mutation.from_string("M198F")
score = scorer.score_mutation(mutation)
print(f"FoldX ΔΔG for {mutation.name}: {score:.3f} kcal/mol")

# Score using mutation name string
score = scorer.score_sequence("M198F")
print(f"ΔΔG: {score:.3f} kcal/mol")
```

## Understanding FoldX Scores

### ΔΔG Interpretation
- **Negative values**: Stabilizing mutation (favorable)
- **Positive values**: Destabilizing mutation (unfavorable)
- **Zero**: Neutral effect on stability

### Typical Ranges
- **< -1.0 kcal/mol**: Strongly stabilizing
- **-1.0 to -0.5 kcal/mol**: Moderately stabilizing
- **-0.5 to 0.5 kcal/mol**: Neutral
- **0.5 to 1.0 kcal/mol**: Moderately destabilizing
- **> 1.0 kcal/mol**: Strongly destabilizing

### Experimental Correlation
- FoldX predictions correlate with experimental ΔΔG values (R ≈ 0.6-0.8)
- Best for single point mutations
- Less accurate for large conformational changes

## Complete Examples

### Example 1: Single Mutation Analysis
```python
#!/usr/bin/env python3
"""
Single mutation analysis with FoldX
"""

import os
import sys
sys.path.insert(0, '/data_test/home/lzzheng/apps/zMutScan/scripts')

from foldx_scorer import FoldXScorer
from base_scorer import Mutation
from utils import pdb2fasta

def analyze_mutation_foldx(pdb_file, mutation_str, output_dir, num_runs=5, is_homodimer=False):
    """Analyze a single mutation with FoldX."""
    
    # Extract sequence from PDB
    sequence = pdb2fasta(pdb_file)
    
    # Initialize scorer
    scorer = FoldXScorer(
        config={
            "args": {
                "num_runs": num_runs,
                "rotabase_location": ""  # Auto-detect
            }
        },
        output_dir=output_dir,
        wt_sequence=sequence,
        pdb_file=pdb_file,
        is_homodimer=is_homodimer
    )
    
    # Create mutation object
    mutation = Mutation.from_string(mutation_str)
    
    # Check if already calculated
    if scorer.result_exists(mutation):
        print(f"Result exists for {mutation_str}, loading...")
        score = scorer.load_result(mutation)
    else:
        # Calculate score
        score = scorer.score_mutation(mutation)
    
    # Interpret result
    interpretation = interpret_foldx_score(score)
    
    return {
        "mutation": mutation_str,
        "ddG": score,
        "interpretation": interpretation,
        "num_runs": num_runs,
        "is_homodimer": is_homodimer
    }

def interpret_foldx_score(ddG):
    """Interpret FoldX ΔΔG score."""
    if ddG < -1.0:
        return "Strongly stabilizing"
    elif ddG < -0.5:
        return "Stabilizing"
    elif -0.5 <= ddG <= 0.5:
        return "Neutral"
    elif ddG <= 1.0:
        return "Destabilizing"
    else:
        return "Strongly destabilizing"

if __name__ == "__main__":
    # Example usage
    result = analyze_mutation_foldx(
        pdb_file="protein.pdb",
        mutation_str="M198F",
        output_dir="./foldx_results",
        num_runs=5,
        is_homodimer=False
    )
    
    print(f"Mutation: {result['mutation']}")
    print(f"FoldX ΔΔG: {result['ddG']:.3f} kcal/mol")
    print(f"Interpretation: {result['interpretation']}")
    print(f"Number of runs: {result['num_runs']}")
    print(f"Homodimer mode: {result['is_homodimer']}")
```

### Example 2: Multiple Mutation Scanning
```python
#!/usr/bin/env python3
"""
Scan multiple mutations with FoldX
"""

import pandas as pd
from tqdm import tqdm

def scan_mutations_foldx(pdb_file, mutations, output_dir, num_runs=5, is_homodimer=False):
    """Scan multiple mutations with FoldX."""
    
    from foldx_scorer import FoldXScorer
    from base_scorer import Mutation
    from utils import pdb2fasta
    
    # Extract sequence
    sequence = pdb2fasta(pdb_file)
    
    # Initialize scorer
    scorer = FoldXScorer(
        config={"args": {"num_runs": num_runs}},
        output_dir=output_dir,
        wt_sequence=sequence,
        pdb_file=pdb_file,
        is_homodimer=is_homodimer
    )
    
    results = []
    
    for mut_str in tqdm(mutations, desc="Scanning mutations with FoldX"):
        try:
            mutation = Mutation.from_string(mut_str)
            
            # Check cache
            if scorer.result_exists(mutation):
                score = scorer.load_result(mutation)
                cached = True
            else:
                score = scorer.score_mutation(mutation)
                cached = False
            
            results.append({
                "mutation": mut_str,
                "ddG": score,
                "cached": cached,
                "error": None
            })
            
        except Exception as e:
            results.append({
                "mutation": mut_str,
                "ddG": None,
                "cached": False,
                "error": str(e)
            })
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Sort by ΔΔG (lower is better for stability)
    df = df.sort_values("ddG", ascending=True)
    
    return df

# Example usage
if __name__ == "__main__":
    # Define mutations to scan
    mutations = [
        "M198F", "M198L", "M198W", "M198Y",
        "K201E", "K201D", "K201R", "K201H",
        "D205N", "D205E", "D205Q", "D205H"
    ]
    
    # Run scan
    results_df = scan_mutations_foldx(
        pdb_file="protein.pdb",
        mutations=mutations,
        output_dir="./foldx_scan_results",
        num_runs=5,
        is_homodimer=False
    )
    
    # Save results
    results_df.to_csv("foldx_scan_results.csv", index=False)
    
    # Display results
    print("\n=== FoldX Mutation Scan Results ===")
    print(f"\nMost stabilizing mutations:")
    print(results_df.head(10))
    
    print(f"\nMost destabilizing mutations:")
    print(results_df.tail(10))
    
    # Statistics
    print(f"\nStatistics:")
    print(f"Total mutations: {len(results_df)}")
    print(f"Stabilizing (ΔΔG < -0.5): {len(results_df[results_df['ddG'] < -0.5])}")
    print(f"Neutral (-0.5 ≤ ΔΔG ≤ 0.5): {len(results_df[(results_df['ddG'] >= -0.5) & (results_df['ddG'] <= 0.5)])}")
    print(f"Destabilizing (ΔΔG > 0.5): {len(results_df[results_df['ddG'] > 0.5])}")
```

### Example 3: Homodimer Analysis
```python
#!/usr/bin/env python3
"""
FoldX analysis for homodimers
"""

def analyze_homodimer_mutations(pdb_file, mutations, output_dir, num_runs=5):
    """Analyze mutations in homodimer context."""
    
    from foldx_scorer import FoldXScorer
    from base_scorer import Mutation
    from utils import pdb2fasta
    
    # Extract sequence (assuming identical chains)
    sequence = pdb2fasta(pdb_file)
    
    # Initialize scorer for homodimer
    scorer = FoldXScorer(
        config={"args": {"num_runs": num_runs}},
        output_dir=output_dir,
        wt_sequence=sequence,
        pdb_file=pdb_file,
        is_homodimer=True
    )
    
    results = []
    
    for mut_str in mutations:
        try:
            mutation = Mutation.from_string(mut_str)
            score = scorer.score_mutation(mutation)
            
            results.append({
                "mutation": mut_str,
                "ddG": score,
                "interpretation": "Homodimer ΔΔG",
                "error": None
            })
            
            print(f"Homodimer mutation {mut_str}: ΔΔG = {score:.3f} kcal/mol")
            
        except Exception as e:
            results.append({
                "mutation": mut_str,
                "ddG": None,
                "interpretation": None,
                "error": str(e)
            })
            print(f"Error for {mut_str}: {e}")
    
    return results

# Compare monomer vs homodimer
def compare_monomer_homodimer(pdb_file, mutation_str, output_base="./comparison"):
    """Compare ΔΔG for monomer vs homodimer."""
    
    from foldx_scorer import FoldXScorer
    from base_scorer import Mutation
    from utils import pdb2fasta
    
    sequence = pdb2fasta(pdb_file)
    mutation = Mutation.from_string(mutation_str)
    
    # Monomer analysis
    monomer_scorer = FoldXScorer(
        config={"args": {"num_runs": 5}},
        output_dir=f"{output_base}/monomer",
        wt_sequence=sequence,
        pdb_file=pdb_file,
        is_homodimer=False
    )
    monomer_ddG = monomer_scorer.score_mutation(mutation)
    
    # Homodimer analysis
    dimer_scorer = FoldXScorer(
        config={"args": {"num_runs": 5}},
        output_dir=f"{output_base}/homodimer",
        wt_sequence=sequence,
        pdb_file=pdb_file,
        is_homodimer=True
    )
    dimer_ddG = dimer_scorer.score_mutation(mutation)
    
    return {
        "mutation": mutation_str,
        "monomer_ddG": monomer_ddG,
        "homodimer_ddG": dimer_ddG,
        "difference": dimer_ddG - monomer_ddG
    }
```

## Advanced Configuration

### Custom Rotabase Location
```python
# Specify custom rotabase file location
scorer = FoldXScorer(
    config={
        "args": {
            "num_runs": 5,
            "rotabase_location": "/path/to/custom/rotabase.txt"
        }
    },
    # ... other parameters
)
```

### Multiple Independent Runs
```python
# Increase number of runs for better statistics
scorer = FoldXScorer(
    config={
        "args": {
            "num_runs": 10,  # More runs = more reliable average
        }
    },
    # ... other parameters
)

# Get run statistics
def get_foldx_statistics(output_dir):
    """Extract statistics from FoldX output."""
    import glob
    import pandas as pd
    
    # Find all output files
    output_files = glob.glob(f"{output_dir}/run_*/Dif_input.fxout")
    
    ddg_values = []
    for file in output_files:
        try:
            with open(file, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith("total"):
                        parts = line.split()
                        if len(parts) >= 3:
                            ddg = float(parts[2])
                            ddg_values.append(ddg)
        except:
            continue
    
    if ddg_values:
        stats = {
            "mean": sum(ddg_values) / len(ddg_values),
            "std": (sum((x - sum(ddg_values)/len(ddg_values))**2 for x in ddg_values) / len(ddg_values))**0.5,
            "min": min(ddg_values),
            "max": max(ddg_values),
            "n": len(ddg_values)
        }
        return stats
    return None
```

## PDB Preparation for FoldX

### Required PDB Format
```python
#!/usr/bin/env python3
"""
Prepare PDB file for FoldX analysis
"""

def prepare_pdb_for_foldx(input_pdb, output_pdb):
    """Clean and prepare PDB file for FoldX."""
    
    from Bio.PDB import PDBParser, PDBIO, Select
    
    class FoldXSelect(Select):
        """Select only standard amino acids for FoldX."""
        def accept_residue(self, residue):
            # Keep only standard amino acids (remove HETATM, water, etc.)
            return residue.id[0] == " "
    
    # Parse input structure
    parser = PDBParser()
    structure = parser.get_structure("input", input_pdb)
    
    # Clean structure
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, FoldXSelect())
    
    print(f"Cleaned PDB saved to: {output_pdb}")
    return output_pdb

def check_pdb_issues(pdb_file):
    """Check for common PDB issues that affect FoldX."""
    
    issues = []
    
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    # Check for missing atoms
    residue_atoms = {}
    for line in lines:
        if line.startswith("ATOM"):
            chain = line[21]
            resnum = int(line[22:26])
            atom = line[12:16].strip()
            
            key = (chain, resnum)
            if key not in residue_atoms:
                residue_atoms[key] = set()
            residue_atoms[key].add(atom)
    
    # Check each residue for backbone atoms
    for (chain, resnum), atoms in residue_atoms.items():
        missing = []
        for required in ["N", "CA", "C", "O"]:
            if required not in atoms:
                missing.append(required)
        
        if missing:
            issues.append(f"Chain {chain} Residue {resnum}: missing {', '.join(missing)}")
    
    return issues
```

## Integration with Other Tools

### Combined Analysis with ESM-IF
```python
#!/usr/bin/env python3
"""
Combine FoldX and ESM-IF scores
"""

def combined_stability_analysis(pdb_file, mutation_str, output_base="./combined"):
    """Run both FoldX and ESM-IF analysis for stability prediction."""
    
    from foldx_scorer import FoldXScorer
    from esmif_scorer import ESMIFScorer
    from base_scorer import Mutation
    from utils import pdb2fasta
    
    sequence = pdb2fasta(pdb_file)
    mutation = Mutation.from_string(mutation_str)
    
    # FoldX analysis (ΔΔG)
    foldx_scorer = FoldXScorer(
        config={"args": {"num_runs": 5}},
        output_dir=f"{output_base}/foldx",
        wt_sequence=sequence,
        pdb_file=pdb_file
    )
    foldx_score = foldx_scorer.score_mutation(mutation)
    
    # ESM-IF analysis (structure compatibility)
    esmif_scorer = ESMIFScorer(
        config={"args": {"model": "esm_if1_gvp4_t16_142M_UR50"}},
        output_dir=f"{output_base}/esmif",
        wt_sequence=sequence,
        pdb_file=pdb_file
    )
    esmif_score = esmif_scorer.score_mutation(mutation)
    
    # Normalize scores (both lower = better)
    # FoldX: negative ΔΔG = stabilizing
    # ESM-IF: negative ΔLoss = better structure compatibility
    
    # Combined score (weighted average)
    combined_score = (foldx_score + esmif_score) / 2
    
    return {
        "mutation": mutation_str,
        "foldx_ddG": foldx_score,
        "esmif_dloss": esmif_score,
        "combined_score": combined_score,
        "interpretation": interpret_combined_score(foldx_score, esmif_score)
    }

def interpret_combined_score(foldx_ddG, esmif_dloss):
    """Interpret combined FoldX and ESM-IF scores."""
    
    # Both negative = favorable
    if foldx_ddG < -0.5 and esmif_dloss < -0.1:
        return "Strongly favorable (stabilizing + good structure fit)"
    elif foldx_ddG < 0 and esmif_dloss < 0:
        return "Favorable"
    elif foldx_ddG > 0.5 and esmif_dloss > 0.1:
        return "Strongly unfavorable (destabilizing + poor structure fit)"
    elif foldx_ddG > 0 and esmif_dloss > 0:
        return "Unfavorable"
    else:
        return "Mixed/neutral"
```

## Troubleshooting

### Common Issues

1. **FoldX not found**
```bash
# Check FoldX installation
which FoldX
which foldx

# Set FoldX path explicitly
export FOLDX_PATH=/path/to/foldx/bin/FoldX
```

2. **Rotabase not found**
```bash
# Find rotabase file
find / -name "rotabase.txt" 2>/dev/null | head -5

# Set rotabase location
export ROTABASE_PATH=/path/to/rotabase.txt
```

3. **PDB format issues**
```python
# Clean PDB file
from Bio.PDB import PDBParser, PDBIO
parser = PDBParser()
structure = parser.get_structure("clean", "protein.pdb")
io = PDBIO()
io.set_structure(structure)
io.save("protein_clean.pdb")
```

4. **Memory errors**
```bash
# Reduce number of runs
python foldx_scorer.py --num_runs 3 protein.pdb M198F ./results

# Clean temporary files
rm -rf ./foldx_results/temp_*
```

### Debug Script
```bash
#!/bin/bash
# debug_foldx.sh

echo "=== FoldX Debug Information ==="

echo -e "\n1. FoldX installation:"
which FoldX 2>/dev/null || which foldx 2>/dev/null || echo "FoldX not found in PATH"

echo -e "\n2. Python environment:"
/data_test/home/lzzheng/.conda/envs/sfct/bin/python --version

echo -e "\n3. FoldX scorer test:"
/data_test/home/lzzheng/.conda/envs/sfct/bin/python -c "
import sys
sys.path.insert(0, '/data_test/home/lzzheng/apps/zMutScan/scripts')
try:
    from foldx_scorer import FoldXScorer
    print('FoldX scorer import successful')
    
    # Test initialization
    scorer = FoldXScorer(
        config={'args': {'num_runs': 1}},
        output_dir='./test_foldx',
        wt_sequence='MRSLTAKEYVEAFKSFLDHSTEHQCMEAFNKQEMPHIMAGLGNGKSTLNVLGVGSGTGEQDLKMIQILQAEHPGVFINAEIIEPNPQHVAAYKELVNQAPGLQNVSFIWHQLTSSEYEQQMKEKSTHKKFDFIHMIQMLYRVEDIPNTIKFFHSCLDHHGKLLIIILSESSGWATLWKKYRHCLPLTDSGHYITSNDIEDILKRIGVEYHVFEFPSGWDITECFIEGDLVGGHMMDFLTGTKNFLGTAPVDLKKRLQEALCQPECSSRKDGRIIFNNNLSMIVVES',
        pdb_file='protein.pdb'
    )
    print('FoldX scorer initialized successfully')
except Exception as e:
    print(f'Error: {e}')
"

echo -e "\n4. PDB file check:"
if [ -f "protein.pdb" ]; then
    echo "PDB file exists"
    head -5 protein.pdb
else
    echo "PDB file not found"
fi
```

## Best Practices

1. **PDB preparation**: Always clean PDB files before analysis
2. **Multiple runs**: Use 5+ runs for reliable ΔΔG estimates
3. **Homodimers**: Use homodimer mode for symmetric complexes
4. **Validation**: Compare with experimental data when available
5. **Combination**: Use FoldX with sequence-based methods for consensus

## References

- **FoldX Paper**: Schymkowitz et al., "The FoldX web server: an online force field", Nucleic Acids Research 2005
- **FoldX Website**: http://foldxsuite.org/
- **ΔΔG Prediction**: Guerois et al., "Predicting changes in the stability of proteins and protein complexes", JMB 2002
- **zMutScan Framework**: Local mutation scanning implementation