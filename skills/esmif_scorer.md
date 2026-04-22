# ESM-IF Scorer Skill

## Overview
ESM-IF (Evolutionary Scale Modeling - Inverse Folding) is a protein language model from Meta AI designed for inverse folding tasks - predicting sequences that fold into a given structure. This skill provides comprehensive guidance on using ESM-IF for mutation scoring through the zMutScan framework.

## Key Features

- **Structure-conditioned**: Uses 3D coordinates to predict sequences
- **Inverse folding**: Optimized for sequence-structure compatibility
- **Single-chain support**: Focuses on individual protein chains
- **Loss-based scoring**: Lower loss = better sequence-structure match

## Installation

### Required Packages
```bash
# Activate sfct environment
source /data_test/home/lzzheng/.conda/envs/sfct/bin/activate

# Install ESM package
pip install fair-esm

# Additional dependencies (usually included with ESM)
pip install torch
pip install biopython
pip install numpy
```

### Verify Installation
```bash
# Test ESM-IF import
/data_test/home/lzzheng/.conda/envs/sfct/bin/python -c "
import esm
import esm.inverse_folding
print('ESM-IF modules available')
print('ESM version:', esm.__version__)
"
```

## Available Models

### Primary Model
- **esm_if1_gvp4_t16_142M_UR50**: 142M parameter model with GVP (Geometric Vector Perceptron) architecture

### Model Characteristics
- **Parameters**: 142 million
- **Architecture**: Transformer with GVP layers
- **Training data**: UR50 (UniRef50)
- **Input**: 3D coordinates (Cα, C, N, O atoms)
- **Output**: Sequence likelihood

## Basic Usage

### Command Line Interface
```bash
# Basic usage
/data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/esmif_scorer.py \
    <pdb_file> \
    <mutation> \
    <output_dir>

# Example
/data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/esmif_scorer.py \
    protein.pdb \
    M198F \
    ./esmif_results

# Specify chain (default is A)
/data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/esmif_scorer.py \
    --chain_id B \
    protein.pdb \
    M198F \
    ./esmif_results
```

### Python API
```python
#!/usr/bin/env python3
from zMutScan.scripts.esmif_scorer import ESMIFScorer
from zMutScan.scripts.base_scorer import Mutation

# Initialize scorer
scorer = ESMIFScorer(
    config={
        "args": {
            "model": "esm_if1_gvp4_t16_142M_UR50",
            "batch_size": 1  # Single mutation mode
        }
    },
    output_dir="./esmif_results",
    wt_sequence="MRSLTAKEYVEAFKSFLDHSTEHQCMEAFNKQEMPHIMAGLGNGKSTLNVLGVGSGTGEQDLKMIQILQAEHPGVFINAEIIEPNPQHVAAYKELVNQAPGLQNVSFIWHQLTSSEYEQQMKEKSTHKKFDFIHMIQMLYRVEDIPNTIKFFHSCLDHHGKLLIIILSESSGWATLWKKYRHCLPLTDSGHYITSNDIEDILKRIGVEYHVFEFPSGWDITECFIEGDLVGGHMMDFLTGTKNFLGTAPVDLKKRLQEALCQPECSSRKDGRIIFNNNLSMIVVES",
    pdb_file="protein.pdb",
    chain_id="A"  # Target chain
)

# Score single mutation
mutation = Mutation.from_string("M198F")
score = scorer.score_mutation(mutation)
print(f"ESM-IF score for {mutation.name}: {score:.6f}")

# Score using mutation name string
score = scorer.score_sequence("M198F")
print(f"Score: {score:.6f}")
```

## How ESM-IF Works

### Input Processing
1. **Structure loading**: PDB file → atomic coordinates
2. **Chain selection**: Extract target chain coordinates
3. **Coordinate extraction**: Cα, C, N, O atoms for each residue
4. **Sequence encoding**: Tokenization of amino acids

### Scoring Process
1. **Wild-type calculation**: Compute loss for native sequence
2. **Mutant calculation**: Compute loss for mutant sequence
3. **Score normalization**: ΔLoss = MutantLoss - WildTypeLoss
4. **Interpretation**: Negative ΔLoss = favorable mutation

## Complete Examples

### Example 1: Single Mutation Analysis
```python
#!/usr/bin/env python3
"""
Single mutation analysis with ESM-IF
"""

import os
import sys
sys.path.insert(0, '/data_test/home/lzzheng/apps/zMutScan/scripts')

from esmif_scorer import ESMIFScorer
from base_scorer import Mutation
from utils import pdb2fasta

def analyze_mutation_esmif(pdb_file, mutation_str, output_dir, chain_id="A"):
    """Analyze a single mutation with ESM-IF."""
    
    # Extract sequence from PDB
    sequence = pdb2fasta(pdb_file)
    
    # Initialize scorer
    scorer = ESMIFScorer(
        config={
            "args": {
                "model": "esm_if1_gvp4_t16_142M_UR50",
                "batch_size": 1
            }
        },
        output_dir=output_dir,
        wt_sequence=sequence,
        pdb_file=pdb_file,
        chain_id=chain_id
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
    interpretation = interpret_esmif_score(score)
    
    return {
        "mutation": mutation_str,
        "score": score,
        "interpretation": interpretation,
        "chain": chain_id,
        "model": "esm_if1_gvp4_t16_142M_UR50"
    }

def interpret_esmif_score(score):
    """Interpret ESM-IF score (ΔLoss)."""
    # Negative score = mutant has lower loss = more compatible with structure
    if score < -0.5:
        return "Strongly favorable (better structure compatibility)"
    elif score < -0.1:
        return "Favorable"
    elif -0.1 <= score <= 0.1:
        return "Neutral"
    elif score <= 0.5:
        return "Unfavorable"
    else:
        return "Strongly unfavorable (worse structure compatibility)"

if __name__ == "__main__":
    # Example usage
    result = analyze_mutation_esmif(
        pdb_file="protein.pdb",
        mutation_str="M198F",
        output_dir="./esmif_results",
        chain_id="A"
    )
    
    print(f"Mutation: {result['mutation']}")
    print(f"Chain: {result['chain']}")
    print(f"ESM-IF ΔLoss: {result['score']:.6f}")
    print(f"Interpretation: {result['interpretation']}")
    print(f"Model: {result['model']}")
```

### Example 2: Multi-Chain Analysis
```python
#!/usr/bin/env python3
"""
Analyze mutations across multiple chains
"""

def analyze_multichain_mutations(pdb_file, mutations_by_chain, output_base="./results"):
    """Analyze mutations in different chains."""
    
    from esmif_scorer import ESMIFScorer
    from base_scorer import Mutation
    from utils import pdb2fasta
    
    # Extract full sequence
    sequence = pdb2fasta(pdb_file)
    
    results = []
    
    for chain_id, mutations in mutations_by_chain.items():
        print(f"\nAnalyzing chain {chain_id}...")
        
        # Initialize scorer for this chain
        scorer = ESMIFScorer(
            config={"args": {"model": "esm_if1_gvp4_t16_142M_UR50"}},
            output_dir=f"{output_base}/chain_{chain_id}",
            wt_sequence=sequence,
            pdb_file=pdb_file,
            chain_id=chain_id
        )
        
        for mut_str in mutations:
            try:
                mutation = Mutation.from_string(mut_str)
                score = scorer.score_mutation(mutation)
                
                results.append({
                    "chain": chain_id,
                    "mutation": mut_str,
                    "score": score,
                    "error": None
                })
                
                print(f"  {mut_str}: {score:.6f}")
                
            except Exception as e:
                results.append({
                    "chain": chain_id,
                    "mutation": mut_str,
                    "score": None,
                    "error": str(e)
                })
                print(f"  {mut_str}: ERROR - {e}")
    
    return results

# Example usage
if __name__ == "__main__":
    # Define mutations by chain
    mutations_by_chain = {
        "A": ["M198F", "K201E", "D205N"],
        "B": ["L302V", "R305K", "E308D"],
        "C": ["S405A", "T408S", "Y411F"]
    }
    
    # Run analysis
    results = analyze_multichain_mutations(
        pdb_file="multichain.pdb",
        mutations_by_chain=mutations_by_chain,
        output_base="./multichain_results"
    )
    
    # Print summary
    print("\n=== Summary ===")
    for result in results:
        if result["error"] is None:
            print(f"Chain {result['chain']} {result['mutation']}: {result['score']:.6f}")
```

### Example 3: Structure Validation Script
```python
#!/usr/bin/env python3
"""
Validate PDB structure for ESM-IF analysis
"""

import esm.inverse_folding
from Bio.PDB import PDBParser

def validate_structure_for_esmif(pdb_file, target_chain="A"):
    """Validate if structure is suitable for ESM-IF."""
    
    issues = []
    
    try:
        # Load structure with ESM-IF utility
        structure = esm.inverse_folding.util.load_structure(pdb_file)
        coords, native_seqs = esm.inverse_folding.multichain_util.extract_coords_from_complex(structure)
        
        # Check if target chain exists
        if target_chain not in coords:
            issues.append(f"Chain {target_chain} not found. Available chains: {list(coords.keys())}")
            return False, issues
        
        # Check coordinate completeness
        chain_coords = coords[target_chain]
        if len(chain_coords) == 0:
            issues.append(f"Chain {target_chain} has no coordinates")
            return False, issues
        
        # Check for missing atoms
        for i, residue_coords in enumerate(chain_coords):
            required_atoms = ["N", "CA", "C"]
            missing_atoms = [atom for atom in required_atoms if atom not in residue_coords]
            if missing_atoms:
                issues.append(f"Residue {i+1} missing atoms: {missing_atoms}")
        
        # Check sequence length
        seq_length = len(native_seqs.get(target_chain, ""))
        coord_length = len(chain_coords)
        if seq_length != coord_length:
            issues.append(f"Sequence length ({seq_length}) != coordinate length ({coord_length})")
        
        return len(issues) == 0, issues
        
    except Exception as e:
        issues.append(f"Structure loading failed: {e}")
        return False, issues

def fix_common_issues(pdb_file, output_file):
    """Attempt to fix common PDB issues for ESM-IF."""
    
    from Bio.PDB import PDBParser, PDBIO
    
    parser = PDBParser()
    structure = parser.get_structure("input", pdb_file)
    
    # Remove hetero atoms and water
    for model in structure:
        for chain in model:
            residues_to_remove = []
            for residue in chain:
                if residue.id[0] != " ":  # Hetero/water
                    residues_to_remove.append(residue.id)
            
            for res_id in residues_to_remove:
                chain.detach_child(res_id)
    
    # Save cleaned structure
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)
    
    return output_file

if __name__ == "__main__":
    pdb_file = "protein.pdb"
    target_chain = "A"
    
    print(f"Validating {pdb_file} for ESM-IF analysis...")
    is_valid, issues = validate_structure_for_esmif(pdb_file, target_chain)
    
    if is_valid:
        print("✓ Structure is valid for ESM-IF")
    else:
        print("✗ Structure has issues:")
        for issue in issues:
            print(f"  - {issue}")
        
        # Try to fix
        print("\nAttempting to fix common issues...")
        fixed_file = fix_common_issues(pdb_file, "protein_cleaned.pdb")
        print(f"Cleaned structure saved to: {fixed_file}")
```

## Advanced Usage

### Batch Processing Multiple Mutations
```python
#!/usr/bin/env python3
"""
Batch process multiple mutations with ESM-IF
"""

import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

def batch_score_esmif(pdb_file, mutations, output_dir, chain_id="A", max_workers=4):
    """Batch score mutations in parallel."""
    
    from esmif_scorer import ESMIFScorer
    from base_scorer import Mutation
    from utils import pdb2fasta
    
    # Extract sequence
    sequence = pdb2fasta(pdb_file)
    
    # Initialize scorer (shared for all workers)
    scorer = ESMIFScorer(
        config={"args": {"model": "esm_if1_gvp4_t16_142M_UR50"}},
        output_dir=output_dir,
        wt_sequence=sequence,
        pdb_file=pdb_file,
        chain_id=chain_id
    )
    
    def score_single(mut_str):
        """Score a single mutation."""
        try:
            mutation = Mutation.from_string(mut_str)
            
            # Check cache first
            if scorer.result_exists(mutation):
                score = scorer.load_result(mutation)
                cached = True
            else:
                score = scorer.score_mutation(mutation)
                cached = False
            
            return {
                "mutation": mut_str,
                "score": score,
                "cached": cached,
                "error": None
            }
        except Exception as e:
            return {
                "mutation": mut_str,
                "score": None,
                "cached": False,
                "error": str(e)
            }
    
    # Process in parallel
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_mut = {executor.submit(score_single, mut): mut for mut in mutations}
        
        for future in as_completed(future_to_mut):
            result = future.result()
            results.append(result)
            print(f"Completed: {result['mutation']} = {result['score']:.6f}")
    
    return pd.DataFrame(results)
```

### Structure-Based Mutation Design
```python
#!/usr/bin/env python3
"""
Structure-based mutation design using ESM-IF
"""

def design_mutations_for_position(pdb_file, position, chain_id="A", output_dir="./design"):
    """Design mutations for a specific position based on structure compatibility."""
    
    import esm
    import esm.inverse_folding
    import torch
    import torch.nn.functional as F
    
    # Load ESM-IF model
    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.eval()
    
    # Load structure
    structure = esm.inverse_folding.util.load_structure(pdb_file)
    coords, native_seqs = esm.inverse_folding.multichain_util.extract_coords_from_complex(structure)
    
    if chain_id not in coords:
        raise ValueError(f"Chain {chain_id} not found")
    
    target_coords = coords[chain_id]
    native_seq = native_seqs.get(chain_id, "")
    
    # Create batch converter
    from esm.inverse_folding.util import CoordBatchConverter
    batch_converter = CoordBatchConverter(alphabet)
    
    # Test all possible mutations at this position
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    results = []
    
    for aa in amino_acids:
        # Skip wild-type amino acid
        if aa == native_seq[position-1]:
            continue
        
        # Create mutant sequence
        mut_seq = native_seq[:position-1] + aa + native_seq[position:]
        
        # Prepare batch
        batch = [(target_coords, None, mut_seq)]
        coords_tensor, confidence, strs, tokens, padding_mask = batch_converter(batch)
        
        # Calculate loss
        with torch.no_grad():
            prev_output_tokens = tokens[:, :-1]
            target = tokens[:, 1:]
            logits, _ = model.forward(coords_tensor, padding_mask, confidence, prev_output_tokens)
            loss = F.cross_entropy(logits, target, reduction='none')
            score = loss.mean().item()
        
        results.append({
            "position": position,
            "mutation": f"{native_seq[position-1]}{position}{aa}",
            "amino_acid": aa,
            "score": score,
            "rank": None
        })
    
    # Sort by score (lower is better)
    results.sort(key=lambda x: x["score"])
    
    # Add ranks
    for i, result in enumerate(results):
        result["rank"] = i + 1
    
    return results
```

## Performance Optimization

### GPU Acceleration
```python
# ESM-IF automatically uses GPU if available
# To check GPU availability:
import torch
print(f"CUDA available: {torch.cuda.is_available()}")
print(f"GPU device: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'None'}")

# Force CPU usage if needed
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""
```

### Memory Management
```bash
# For large structures, reduce batch size
python esmif_scorer.py --batch_size 1 protein.pdb M198F ./results

# Monitor GPU memory
nvidia-smi -l 1  # Monitor every second
```

### Caching Strategy
```python
# ESM-IF scorer implements automatic caching
# Cache location: <output_dir>/results.json

# To manually manage cache:
import json

def load_cache(output_dir):
    cache_file = os.path.join(output_dir, "results.json")
    if os.path.exists(cache_file):
        with open(cache_file, 'r') as f:
            return json.load(f)
    return {}

def save_to_cache(output_dir, mutation, score, metadata):
    cache_file = os.path.join(output_dir, "results.json")
    cache = load_cache(output_dir)
    
    cache[mutation.name] = {
        "score": score,
        "metadata": metadata,
        "timestamp": datetime.now().isoformat()
    }
    
    with open(cache_file, 'w') as f:
        json.dump(cache, f, indent=2)
```

## Troubleshooting

### Common Issues

1. **Structure loading errors**
```bash
# Check PDB format
head -50 protein.pdb

# Try cleaning the PDB
/data_test/home/lzzheng/.conda/envs/sfct/bin/python -c "
from Bio.PDB import PDBParser, PDBIO
parser = PDBParser()
structure = parser.get_structure('clean', 'protein.pdb')
io = PDBIO()
io.set_structure(structure)
io.save('protein_clean.pdb')
print('Cleaned PDB saved')
"
```

2. **Chain not found**
```python
# List available chains
import esm.inverse_folding
structure = esm.inverse_folding.util.load_structure("protein.pdb")
coords, native_seqs = esm.inverse_folding.multichain_util.extract_coords_from_complex(structure)
print(f"Available chains: {list(coords.keys())}")
```

3. **Memory errors**
```bash
# Use CPU instead of GPU
CUDA_VISIBLE_DEVICES="" python esmif_scorer.py protein.pdb M198F ./results

# Reduce model precision
python esmif_scorer.py --precision fp16 protein.pdb M198F ./results
```

### Debug Script
```bash
#!/bin/bash
# debug_esmif.sh

echo "=== ESM-IF Debug Information ==="

echo -e "\n1. Python environment:"
/data_test/home/lzzheng/.conda/envs/sfct/bin/python --version

echo -e "\n2. ESM-IF availability:"
/data_test/home/lzzheng/.conda/envs/sfct/bin/python -c "
import esm
import esm.inverse_folding
print('ESM version:', esm.__version__)
print('Inverse folding module available')

# Test model loading
try:
    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    print('ESM-IF model loaded successfully')
    print(f'Model parameters: {sum(p.numel() for p in model.parameters()):,}')
except Exception as e:
    print(f'Model loading failed: {e}')
"

echo -e "\n3. PDB structure test:"
/data_test/home/lzzheng/.conda/envs/sfct/bin/python -c "
import esm.inverse_folding
try:
    structure = esm.inverse_folding.util.load_structure('protein.pdb')
    coords, seqs = esm.inverse_folding.multichain_util.extract_coords_from_complex(structure)
    print(f'Structure loaded successfully')
    print(f'Chains: {list(coords.keys())}')
    for chain, seq in seqs.items():
        print(f'  Chain {chain}: {len(seq)} residues')
except Exception as e:
    print(f'Structure loading failed: {e}')
"
```

## Best Practices

1. **Structure preparation**: Clean PDB files before analysis
2. **Chain selection**: Specify correct chain ID for multi-chain proteins
3. **Score interpretation**: Negative ΔLoss = favorable mutation
4. **Validation**: Always validate mutations against wild-type sequence
5. **Combination**: Use ESM-IF with other tools (FoldX, VESPA) for consensus

## References

- **ESM-IF Paper**: Hsu et al., "Learning inverse folding from millions of predicted structures", ICML 2022
- **ESM GitHub**: https://github.com/facebookresearch/esm
- **Inverse Folding Tutorial**: https://github.com/facebookresearch/esm/tree/main/examples/inverse_folding
- **zMutScan Framework**: Local mutation scanning implementation