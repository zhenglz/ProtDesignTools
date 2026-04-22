# VESPA Scorer Skill

## Overview
VESPA (Variant Effect Score Prediction without Alignments) is a method that uses protein language model embeddings (ESM2) to predict mutation effects without requiring multiple sequence alignments. This skill provides comprehensive guidance on using the VESPA implementation in the zMutScan framework.

## Key Features

- **No MSA required**: Uses single-sequence embeddings
- **ESM2-based**: Leverages state-of-the-art protein language models
- **Fast prediction**: Suitable for high-throughput mutation scanning
- **Structure-aware**: Can incorporate structural information when available

## Installation

### Required Packages
```bash
# Activate sfct environment
source /data_test/home/lzzheng/.conda/envs/sfct/bin/activate

# Install core dependencies
pip install torch
pip install fair-esm
pip install biopython
pip install pandas
pip install numpy
```

### Verify Installation
```bash
# Test VESPA imports
/data_test/home/lzzheng/.conda/envs/sfct/bin/python -c "
import sys
sys.path.insert(0, '/data_test/home/lzzheng/apps/zMutScan/scripts')
try:
    from vespa_scorer import VESPAScorer
    print('VESPA scorer import successful')
except ImportError as e:
    print(f'Import error: {e}')
"
```

## Basic Usage

### Command Line Interface
```bash
# Basic usage
/data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/vespa_scorer.py \
    <pdb_file> \
    <mutation> \
    <output_dir>

# Example with ESM2 650M model
/data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/vespa_scorer.py \
    protein.pdb \
    M198F \
    ./vespa_results

# Specify model size
/data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/vespa_scorer.py \
    --model esm2_150M \
    protein.pdb \
    M198F \
    ./vespa_results
```

### Python API
```python
#!/usr/bin/env python3
from zMutScan.scripts.vespa_scorer import VESPAScorer
from zMutScan.scripts.base_scorer import Mutation

# Initialize scorer
scorer = VESPAScorer(
    config={
        "args": {
            "model": "esm2_650M",  # or "esm2_35M", "esm2_150M", "esm2_3B"
            "local_model_path": "/path/to/custom/model.pt"  # Optional
        }
    },
    output_dir="./vespa_results",
    wt_sequence="MRSLTAKEYVEAFKSFLDHSTEHQCMEAFNKQEMPHIMAGLGNGKSTLNVLGVGSGTGEQDLKMIQILQAEHPGVFINAEIIEPNPQHVAAYKELVNQAPGLQNVSFIWHQLTSSEYEQQMKEKSTHKKFDFIHMIQMLYRVEDIPNTIKFFHSCLDHHGKLLIIILSESSGWATLWKKYRHCLPLTDSGHYITSNDIEDILKRIGVEYHVFEFPSGWDITECFIEGDLVGGHMMDFLTGTKNFLGTAPVDLKKRLQEALCQPECSSRKDGRIIFNNNLSMIVVES",
    pdb_file="protein.pdb",
    chain_id="A"  # Optional, default is "A"
)

# Score single mutation
mutation = Mutation.from_string("M198F")
score = scorer.score_mutation(mutation)
print(f"VESPA score for {mutation.name}: {score:.6f}")

# Score using mutation name string
score = scorer.score_sequence("M198F")
print(f"Score: {score:.6f}")
```

## Advanced Configuration

### Model Selection
```python
# Different ESM2 models available
MODEL_OPTIONS = {
    "esm2_35M": "esm2_t12_35M_UR50D",      # 35 million parameters
    "esm2_150M": "esm2_t30_150M_UR50D",    # 150 million parameters
    "esm2_650M": "esm2_t33_650M_UR50D",    # 650 million parameters (recommended)
    "esm2_3B": "esm2_t36_3B_UR50D",        # 3 billion parameters
}

# Initialize with specific model
scorer = VESPAScorer(
    config={
        "args": {
            "model": "esm2_3B",  # Use largest model
            "batch_size": 4,     # Batch size for processing
        }
    },
    # ... other parameters
)
```

### Local Model Files
```python
# Use pre-downloaded model files
scorer = VESPAScorer(
    config={
        "args": {
            "model": "esm2_650M",
            "local_model_path": "~/.cache/torch/hub/checkpoints/esm2_t33_650M_UR50D.pt"
        }
    },
    # ... other parameters
)
```

## Complete Examples

### Example 1: Single Mutation Analysis
```python
#!/usr/bin/env python3
"""
Single mutation analysis with VESPA
"""

import os
import sys
sys.path.insert(0, '/data_test/home/lzzheng/apps/zMutScan/scripts')

from vespa_scorer import VESPAScorer
from base_scorer import Mutation
from utils import pdb2fasta

def analyze_mutation(pdb_file, mutation_str, output_dir, model="esm2_650M"):
    """Analyze a single mutation."""
    
    # Extract sequence from PDB
    sequence = pdb2fasta(pdb_file)
    
    # Initialize scorer
    scorer = VESPAScorer(
        config={"args": {"model": model}},
        output_dir=output_dir,
        wt_sequence=sequence,
        pdb_file=pdb_file,
        chain_id="A"
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
    interpretation = interpret_vespa_score(score)
    
    return {
        "mutation": mutation_str,
        "score": score,
        "interpretation": interpretation,
        "model": model
    }

def interpret_vespa_score(score):
    """Interpret VESPA score."""
    if score < -2:
        return "Strongly favorable"
    elif score < 0:
        return "Favorable"
    elif score == 0:
        return "Neutral"
    elif score < 2:
        return "Unfavorable"
    else:
        return "Strongly unfavorable"

if __name__ == "__main__":
    # Example usage
    result = analyze_mutation(
        pdb_file="protein.pdb",
        mutation_str="M198F",
        output_dir="./results",
        model="esm2_650M"
    )
    
    print(f"Mutation: {result['mutation']}")
    print(f"VESPA score: {result['score']:.6f}")
    print(f"Interpretation: {result['interpretation']}")
    print(f"Model used: {result['model']}")
```

### Example 2: Multiple Mutation Scanning
```python
#!/usr/bin/env python3
"""
Scan multiple mutations with VESPA
"""

import pandas as pd
from tqdm import tqdm

def scan_mutations(pdb_file, mutations, output_dir, model="esm2_650M"):
    """Scan multiple mutations."""
    
    from vespa_scorer import VESPAScorer
    from base_scorer import Mutation
    from utils import pdb2fasta
    
    # Extract sequence
    sequence = pdb2fasta(pdb_file)
    
    # Initialize scorer
    scorer = VESPAScorer(
        config={"args": {"model": model}},
        output_dir=output_dir,
        wt_sequence=sequence,
        pdb_file=pdb_file
    )
    
    results = []
    
    for mut_str in tqdm(mutations, desc="Scanning mutations"):
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
                "score": score,
                "cached": cached,
                "error": None
            })
            
        except Exception as e:
            results.append({
                "mutation": mut_str,
                "score": None,
                "cached": False,
                "error": str(e)
            })
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Sort by score (lower is better for VESPA)
    df = df.sort_values("score", ascending=True)
    
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
    results_df = scan_mutations(
        pdb_file="protein.pdb",
        mutations=mutations,
        output_dir="./scan_results",
        model="esm2_650M"
    )
    
    # Save results
    results_df.to_csv("vespa_scan_results.csv", index=False)
    
    # Display top mutations
    print("\nTop favorable mutations:")
    print(results_df.head(10))
    
    print("\nTop unfavorable mutations:")
    print(results_df.tail(10))
```

### Example 3: Command Line Wrapper Script
```bash
#!/bin/bash
# vespa_scan.sh - Batch VESPA scanning wrapper

set -e

PDB_FILE=$1
MUTATION_LIST=$2
OUTPUT_DIR=$3
MODEL=${4:-"esm2_650M"}

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Extract sequence
SEQUENCE=$(/data_test/home/lzzheng/.conda/envs/sfct/bin/python -c "
import sys
sys.path.insert(0, '/data_test/home/lzzheng/apps/zMutScan/scripts')
from utils import pdb2fasta
print(pdb2fasta('$PDB_FILE'))
")

echo "Protein sequence length: ${#SEQUENCE}"
echo "Model: $MODEL"
echo "Output directory: $OUTPUT_DIR"

# Process each mutation
while read -r MUTATION; do
    if [[ -z "$MUTATION" ]] || [[ "$MUTATION" == \#* ]]; then
        continue
    fi
    
    echo "Processing mutation: $MUTATION"
    
    /data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/vespa_scorer.py \
        "$PDB_FILE" \
        "$MUTATION" \
        "$OUTPUT_DIR" \
        2>&1 | tee "$OUTPUT_DIR/${MUTATION}.log"
    
    # Check exit status
    if [ $? -eq 0 ]; then
        echo "✓ Completed: $MUTATION"
    else
        echo "✗ Failed: $MUTATION"
    fi
    
done < "$MUTATION_LIST"

echo "VESPA scanning complete!"
```

## Performance Optimization

### GPU Acceleration
```python
# VESPA automatically uses GPU if available
# To force CPU usage:
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""

# Or specify device in scorer initialization
scorer = VESPAScorer(
    config={
        "args": {
            "model": "esm2_650M",
            "device": "cpu"  # Force CPU
        }
    },
    # ... other parameters
)
```

### Batch Processing
```python
# VESPA supports batch processing for multiple mutations
scorer = VESPAScorer(
    config={
        "args": {
            "model": "esm2_650M",
            "batch_size": 8  # Process 8 mutations simultaneously
        }
    },
    # ... other parameters
)
```

### Caching Results
```python
# VESPA automatically caches results in output_dir
# To clear cache:
import shutil
shutil.rmtree("./vespa_results")

# Or disable caching:
scorer = VESPAScorer(
    config={
        "args": {
            "model": "esm2_650M",
            "cache_results": False  # Disable caching
        }
    },
    output_dir="./temp_results",  # Still need output_dir for temp files
    # ... other parameters
)
```

## Integration with Other Tools

### Combined Analysis with FoldX
```python
#!/usr/bin/env python3
"""
Combine VESPA and FoldX scores
"""

from vespa_scorer import VESPAScorer
from foldx_scorer import FoldXScorer
from base_scorer import Mutation

def combined_analysis(pdb_file, mutation_str, output_base="./results"):
    """Run both VESPA and FoldX analysis."""
    
    from utils import pdb2fasta
    sequence = pdb2fasta(pdb_file)
    mutation = Mutation.from_string(mutation_str)
    
    # VESPA analysis
    vespa_scorer = VESPAScorer(
        config={"args": {"model": "esm2_650M"}},
        output_dir=f"{output_base}/vespa",
        wt_sequence=sequence,
        pdb_file=pdb_file
    )
    vespa_score = vespa_scorer.score_mutation(mutation)
    
    # FoldX analysis
    foldx_scorer = FoldXScorer(
        config={"args": {"num_runs": 5}},
        output_dir=f"{output_base}/foldx",
        wt_sequence=sequence,
        pdb_file=pdb_file
    )
    foldx_score = foldx_scorer.score_mutation(mutation)
    
    # Combine scores (normalized)
    combined_score = (vespa_score + foldx_score) / 2
    
    return {
        "mutation": mutation_str,
        "vespa_score": vespa_score,
        "foldx_score": foldx_score,
        "combined_score": combined_score
    }
```

### Visualization
```python
import matplotlib.pyplot as plt
import seaborn as sns

def visualize_vespa_scores(results_df):
    """Visualize VESPA scores."""
    
    plt.figure(figsize=(10, 6))
    
    # Sort by score
    df_sorted = results_df.sort_values("score")
    
    # Create bar plot
    plt.bar(range(len(df_sorted)), df_sorted["score"])
    plt.xticks(range(len(df_sorted)), df_sorted["mutation"], rotation=45, ha='right')
    plt.ylabel("VESPA Score (lower is better)")
    plt.title("VESPA Mutation Scores")
    plt.axhline(y=0, color='r', linestyle='-', alpha=0.3)
    
    # Color code by score
    for i, score in enumerate(df_sorted["score"]):
        if score < -1:
            plt.gca().get_children()[i].set_color('green')
        elif score > 1:
            plt.gca().get_children()[i].set_color('red')
    
    plt.tight_layout()
    plt.savefig("vespa_scores.png", dpi=300)
    plt.show()
```

## Troubleshooting

### Common Issues

1. **Model loading errors**
```bash
# Check if model files exist
ls ~/.cache/torch/hub/checkpoints/

# Download model manually
wget https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt \
     -P ~/.cache/torch/hub/checkpoints/
```

2. **Out of memory**
```bash
# Use smaller model
python vespa_scorer.py --model esm2_150M protein.pdb M198F ./results

# Reduce batch size
export VESPA_BATCH_SIZE=1
```

3. **PDB parsing errors**
```python
# Check PDB file
from Bio.PDB import PDBParser
parser = PDBParser()
structure = parser.get_structure("test", "protein.pdb")
print(f"Chains: {[chain.id for chain in structure.get_chains()]}")
```

### Debug Script
```bash
#!/bin/bash
# debug_vespa.sh

echo "=== VESPA Debug Information ==="

echo -e "\n1. Python environment:"
/data_test/home/lzzheng/.conda/envs/sfct/bin/python --version

echo -e "\n2. Package versions:"
/data_test/home/lzzheng/.conda/envs/sfct/bin/python -c "
import torch, esm, sys
print(f'Torch: {torch.__version__}')
print(f'ESM: {esm.__version__}')
print(f'Python: {sys.version}')
"

echo -e "\n3. Model availability:"
/data_test/home/lzzheng/.conda/envs/sfct/bin/python -c "
import esm
try:
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    print('ESM2 650M model loaded successfully')
except Exception as e:
    print(f'Model loading failed: {e}')
"

echo -e "\n4. VESPA import test:"
/data_test/home/lzzheng/.conda/envs/sfct/bin/python -c "
import sys
sys.path.insert(0, '/data_test/home/lzzheng/apps/zMutScan/scripts')
try:
    from vespa_scorer import VESPAScorer
    print('VESPA scorer import successful')
except Exception as e:
    print(f'Import error: {e}')
"
```

## Best Practices

1. **Model selection**: Use `esm2_650M` for best accuracy/performance balance
2. **Batch processing**: Use batch_size > 1 for multiple mutations
3. **Caching**: Enable caching for reproducible results
4. **Validation**: Always validate mutations against wild-type sequence
5. **Interpretation**: Consider VESPA scores in context with other metrics

## References

- **VESPA Paper**: Marquet et al., "Embeddings from protein language models predict conservation and variant effects", Human Genetics 2021
- **ESM2 Paper**: Lin et al., "Language models of protein sequences at the scale of evolution enable accurate structure prediction", bioRxiv 2022
- **zMutScan Repository**: Local mutation scanning framework
- **ESM GitHub**: https://github.com/facebookresearch/esm