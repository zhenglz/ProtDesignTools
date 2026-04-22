# ESM2 Scorer Skill

## Overview
ESM2 (Evolutionary Scale Modeling 2) is a protein language model from Meta AI that can be used for mutation effect prediction. This skill provides guidance on using ESM2 models for scoring protein mutations through the VESPA framework or direct ESM2 inference.

## Available Models

### ESM2 Model Sizes
1. **ESM2 35M** (`esm2_t12_35M_UR50D`) - Smallest, fastest
2. **ESM2 150M** (`esm2_t30_150M_UR50D`) - Balanced
3. **ESM2 650M** (`esm2_t33_650M_UR50D`) - Recommended for most tasks
4. **ESM2 3B** (`esm2_t36_3B_UR50D`) - Largest, most accurate

### Local Model Paths
```bash
# Default cache locations
~/.cache/torch/hub/checkpoints/esm2_t12_35M_UR50D.pt
~/.cache/torch/hub/checkpoints/esm2_t30_150M_UR50D.pt
~/.cache/torch/hub/checkpoints/esm2_t33_650M_UR50D.pt
~/.cache/torch/hub/checkpoints/esm2_t36_3B_UR50D.pt
```

## Installation

### Python Environment
```bash
# Activate the sfct environment
source /data_test/home/lzzheng/.conda/envs/sfct/bin/activate

# Install required packages
pip install torch
pip install fair-esm  # ESM package
pip install biopython
```

### Check Installation
```bash
# Test ESM2 import
/data_test/home/lzzheng/.conda/envs/sfct/bin/python -c "import esm; print('ESM version:', esm.__version__)"
```

## Usage Methods

### Method 1: Using VESPA Scorer (Recommended)
VESPA (Variant Effect Score Prediction without Alignments) uses ESM2 embeddings to predict mutation effects.

```bash
# Command line usage
/data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/vespa_scorer.py \
    <pdb_file> \
    <mutation> \
    <output_dir>

# Example
/data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/vespa_scorer.py \
    protein.pdb \
    M198F \
    ./vespa_results
```

### Method 2: Direct ESM2 Inference
```python
#!/usr/bin/env python3
import esm
import torch

# Load ESM2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model = model.eval()

# Prepare sequence
sequence = "MRSLTAKEYVEAFKSFLDHSTEHQCMEAFNKQEMPHIMAGLGNGKSTLNVLGVGSGTGEQDLKMIQILQAEHPGVFINAEIIEPNPQHVAAYKELVNQAPGLQNVSFIWHQLTSSEYEQQMKEKSTHKKFDFIHMIQMLYRVEDIPNTIKFFHSCLDHHGKLLIIILSESSGWATLWKKYRHCLPLTDSGHYITSNDIEDILKRIGVEYHVFEFPSGWDITECFIEGDLVGGHMMDFLTGTKNFLGTAPVDLKKRLQEALCQPECSSRKDGRIIFNNNLSMIVVES"
batch_converter = alphabet.get_batch_converter()

# Prepare data
data = [("protein1", sequence)]
batch_labels, batch_strs, batch_tokens = batch_converter(data)

# Get embeddings
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33])
    embeddings = results["representations"][33]
    
print(f"Embedding shape: {embeddings.shape}")
```

## Python Script Examples

### Complete ESM2 Mutation Scorer
```python
#!/usr/bin/env python3
"""
ESM2 Mutation Scorer - Direct ESM2-based mutation scoring
"""

import os
import sys
import argparse
import torch
import numpy as np

def load_esm2_model(model_name="esm2_t33_650M_UR50D"):
    """Load ESM2 model."""
    import esm
    
    model_map = {
        "esm2_35M": esm.pretrained.esm2_t12_35M_UR50D,
        "esm2_150M": esm.pretrained.esm2_t30_150M_UR50D,
        "esm2_650M": esm.pretrained.esm2_t33_650M_UR50D,
        "esm2_3B": esm.pretrained.esm2_t36_3B_UR50D,
    }
    
    if model_name in model_map:
        model, alphabet = model_map[model_name]()
    else:
        # Default to 650M
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    
    model = model.eval()
    return model, alphabet

def calculate_mutation_score(sequence, mutation, model, alphabet):
    """
    Calculate ESM2 score for a mutation.
    
    Args:
        sequence: Wild-type protein sequence
        mutation: Mutation string like "M198F"
        model: ESM2 model
        alphabet: ESM2 alphabet
    
    Returns:
        float: Mutation score (negative log likelihood difference)
    """
    # Parse mutation
    wt_aa = mutation[0]
    position = int(mutation[1:-1])
    mut_aa = mutation[-1]
    
    # Validate
    if position < 1 or position > len(sequence):
        raise ValueError(f"Position {position} out of range (1-{len(sequence)})")
    
    if sequence[position-1] != wt_aa:
        raise ValueError(f"Wild-type mismatch: expected {wt_aa}, found {sequence[position-1]}")
    
    # Create mutant sequence
    mut_sequence = sequence[:position-1] + mut_aa + sequence[position:]
    
    # Prepare batch converter
    batch_converter = alphabet.get_batch_converter()
    
    # Calculate wild-type log likelihood
    wt_data = [("wt", sequence)]
    wt_labels, wt_strs, wt_tokens = batch_converter(wt_data)
    
    with torch.no_grad():
        wt_results = model(wt_tokens, repr_layers=[33])
        wt_logits = wt_results["logits"]
        
        # Get position-specific logits
        pos_logits = wt_logits[0, position]  # Shape: (vocab_size,)
        
        # Convert amino acids to token indices
        wt_token_idx = alphabet.tok_to_idx[wt_aa]
        mut_token_idx = alphabet.tok_to_idx[mut_aa]
        
        # Calculate log likelihoods
        wt_log_prob = torch.log_softmax(pos_logits, dim=0)[wt_token_idx].item()
        mut_log_prob = torch.log_softmax(pos_logits, dim=0)[mut_token_idx].item()
        
        # Score is negative log likelihood difference
        score = -(mut_log_prob - wt_log_prob)
    
    return score

def main():
    parser = argparse.ArgumentParser(description="ESM2 Mutation Scorer")
    parser.add_argument("--sequence", required=True, help="Protein sequence")
    parser.add_argument("--mutation", required=True, help="Mutation like M198F")
    parser.add_argument("--model", default="esm2_650M", 
                       choices=["esm2_35M", "esm2_150M", "esm2_650M", "esm2_3B"],
                       help="ESM2 model to use")
    
    args = parser.parse_args()
    
    # Load model
    print(f"Loading ESM2 model: {args.model}")
    model, alphabet = load_esm2_model(args.model)
    
    # Calculate score
    try:
        score = calculate_mutation_score(args.sequence, args.mutation, model, alphabet)
        print(f"Mutation {args.mutation}: {score:.6f}")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
```

### Batch Mutation Scoring
```python
#!/usr/bin/env python3
"""
Batch ESM2 mutation scoring
"""

import pandas as pd
from tqdm import tqdm

def batch_score_mutations(sequence, mutations, model_name="esm2_650M"):
    """Score multiple mutations."""
    import esm
    import torch
    
    # Load model
    model, alphabet = load_esm2_model(model_name)
    batch_converter = alphabet.get_batch_converter()
    
    results = []
    
    for mutation in tqdm(mutations, desc="Scoring mutations"):
        try:
            score = calculate_mutation_score(sequence, mutation, model, alphabet)
            results.append({
                "mutation": mutation,
                "score": score,
                "error": None
            })
        except Exception as e:
            results.append({
                "mutation": mutation,
                "score": None,
                "error": str(e)
            })
    
    return pd.DataFrame(results)
```

## Integration with zMutScan Framework

### Using the VESPA Scorer Class
```python
from zMutScan.scripts.vespa_scorer import VESPAScorer
from zMutScan.scripts.base_scorer import Mutation

# Initialize scorer
scorer = VESPAScorer(
    config={
        "args": {
            "model": "esm2_650M",
            "local_model_path": "/path/to/local/model.pt"  # Optional
        }
    },
    output_dir="./results",
    wt_sequence="MRSLTAKEYVEAFKSFLDHSTEHQCMEAFNKQEMPHIMAGLGNGKSTLNVLGVGSGTGEQDLKMIQILQAEHPGVFINAEIIEPNPQHVAAYKELVNQAPGLQNVSFIWHQLTSSEYEQQMKEKSTHKKFDFIHMIQMLYRVEDIPNTIKFFHSCLDHHGKLLIIILSESSGWATLWKKYRHCLPLTDSGHYITSNDIEDILKRIGVEYHVFEFPSGWDITECFIEGDLVGGHMMDFLTGTKNFLGTAPVDLKKRLQEALCQPECSSRKDGRIIFNNNLSMIVVES",
    pdb_file="protein.pdb"
)

# Score mutation
mutation = Mutation.from_string("M198F")
score = scorer.score_mutation(mutation)
print(f"VESPA score for {mutation.name}: {score}")
```

### Command Line Wrapper Script
```bash
#!/bin/bash
# esm2_score.sh - Wrapper for ESM2 scoring

set -e

PDB_FILE=$1
MUTATION=$2
OUTPUT_DIR=$3
MODEL=${4:-"esm2_650M"}

# Extract sequence from PDB
SEQUENCE=$(/data_test/home/lzzheng/.conda/envs/sfct/bin/python -c "
import sys
sys.path.insert(0, '/data_test/home/lzzheng/apps/zMutScan/scripts')
from utils import pdb2fasta
print(pdb2fasta('$PDB_FILE'))
")

# Run scoring
/data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/vespa_scorer.py \
    --model "$MODEL" \
    "$PDB_FILE" \
    "$MUTATION" \
    "$OUTPUT_DIR"

echo "ESM2 scoring complete. Results in $OUTPUT_DIR"
```

## Performance Tips

### GPU Acceleration
```python
# Use GPU if available
import torch
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = model.to(device)
```

### Batch Processing
```python
# Process multiple mutations in batch
def batch_score(sequences, mutations, model, alphabet):
    batch_converter = alphabet.get_batch_converter()
    
    # Prepare batch data
    batch_data = []
    for seq, mut in zip(sequences, mutations):
        batch_data.append((f"mut_{mut}", seq))
    
    # Convert batch
    batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
    
    # Move to GPU if available
    if torch.cuda.is_available():
        batch_tokens = batch_tokens.cuda()
    
    # Forward pass
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33])
    
    return results
```

### Memory Management
```bash
# For large models (3B), use CPU or reduce batch size
export CUDA_VISIBLE_DEVICES=""  # Force CPU usage
# or
python script.py --model esm2_150M  # Use smaller model
```

## Interpretation of Scores

### Score Meaning
- **Negative values**: Mutation is favorable (lower is better)
- **Positive values**: Mutation is unfavorable
- **Zero**: Neutral mutation

### Typical Ranges
- **< -2**: Strongly favorable mutation
- **-2 to 0**: Mildly favorable to neutral
- **0 to 2**: Mildly unfavorable
- **> 2**: Strongly unfavorable

## Troubleshooting

### Common Issues

1. **Model loading fails**
```bash
# Check if model files exist
ls ~/.cache/torch/hub/checkpoints/

# Download models manually
wget https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt -P ~/.cache/torch/hub/checkpoints/
```

2. **Out of memory**
```bash
# Use smaller model
python script.py --model esm2_150M

# Use CPU instead of GPU
export CUDA_VISIBLE_DEVICES=""
```

3. **Import errors**
```bash
# Reinstall ESM
pip uninstall fair-esm -y
pip install fair-esm

# Check Python version
python --version  # Should be 3.6+
```

### Debugging Script
```bash
#!/bin/bash
# debug_esm2.sh

echo "=== ESM2 Debug Information ==="
echo "Python version:"
/data_test/home/lzzheng/.conda/envs/sfct/bin/python --version

echo -e "\nESM import test:"
/data_test/home/lzzheng/.conda/envs/sfct/bin/python -c "
try:
    import esm
    print('ESM version:', esm.__version__)
    
    # Test model loading
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    print('Model loaded successfully')
    print(f'Model parameters: {sum(p.numel() for p in model.parameters()):,}')
except Exception as e:
    print(f'Error: {e}')
"

echo -e "\nTorch CUDA available:"
/data_test/home/lzzheng/.conda/envs/sfct/bin/python -c "import torch; print(torch.cuda.is_available())"
```

## References

- **ESM2 Paper**: Lin et al., "Language models of protein sequences at the scale of evolution enable accurate structure prediction", bioRxiv 2022
- **VESPA Paper**: Marquet et al., "Embeddings from protein language models predict conservation and variant effects", Human Genetics 2021
- **ESM GitHub**: https://github.com/facebookresearch/esm
- **Model Downloads**: https://github.com/facebookresearch/esm#available