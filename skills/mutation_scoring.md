# Mutation Scoring Tools Skill

## Overview
This skill provides a comprehensive guide to using various mutation scoring tools available in the zMutScan framework. It covers ESM2, VESPA, ESM-IF, and FoldX for protein mutation effect prediction.

## Quick Comparison

| Tool | Type | Input | Output | Interpretation | Best For |
|------|------|-------|--------|----------------|----------|
| **ESM2** | Language Model | Sequence | Log likelihood | Lower = better | Sequence-only prediction |
| **VESPA** | Embedding-based | Sequence/Structure | ΔLog odds | Lower = better | Fast, no MSA required |
| **ESM-IF** | Inverse Folding | Structure | ΔLoss | Lower = better | Structure compatibility |
| **FoldX** | Force Field | Structure | ΔΔG (kcal/mol) | Lower = better | Stability prediction |

## Installation Summary

### Python Environment
```bash
# Activate sfct environment
source /data_test/home/lzzheng/.conda/envs/sfct/bin/activate

# Common dependencies
pip install torch
pip install fair-esm
pip install biopython
pip install pandas
pip install numpy
```

### Tool-Specific Requirements
```bash
# ESM2/VESPA/ESM-IF
pip install fair-esm  # Meta ESM package

# FoldX
# Requires FoldX binary (http://foldxsuite.org/)
# and rotabase.txt file
```

## Quick Start Examples

### Single Mutation Analysis
```bash
# ESM2 (via VESPA)
/data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/vespa_scorer.py \
    protein.pdb M198F ./results

# ESM-IF
/data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/esmif_scorer.py \
    protein.pdb M198F ./results

# FoldX
/data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/foldx_scorer.py \
    protein.pdb M198F ./results
```

### Batch Processing Script
```bash
#!/bin/bash
# batch_score.sh - Score mutations with multiple tools

set -e

PDB_FILE=$1
MUTATION_LIST=$2
OUTPUT_BASE=$3

# Process each mutation
while read -r MUTATION; do
    if [[ -z "$MUTATION" ]] || [[ "$MUTATION" == \#* ]]; then
        continue
    fi
    
    echo "=== Processing $MUTATION ==="
    
    # VESPA
    echo "Running VESPA..."
    /data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/vespa_scorer.py \
        "$PDB_FILE" "$MUTATION" "$OUTPUT_BASE/vespa" 2>&1 | tail -5
    
    # ESM-IF
    echo "Running ESM-IF..."
    /data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/esmif_scorer.py \
        "$PDB_FILE" "$MUTATION" "$OUTPUT_BASE/esmif" 2>&1 | tail -5
    
    # FoldX
    echo "Running FoldX..."
    /data_test/home/lzzheng/.conda/envs/sfct/bin/python /data_test/home/lzzheng/apps/zMutScan/scripts/foldx_scorer.py \
        "$PDB_FILE" "$MUTATION" "$OUTPUT_BASE/foldx" 2>&1 | tail -5
    
    echo ""
    
done < "$MUTATION_LIST"

echo "Batch scoring complete!"
```

## Python Integration

### Unified Scoring Interface
```python
#!/usr/bin/env python3
"""
Unified mutation scoring interface
"""

import os
import sys
sys.path.insert(0, '/data_test/home/lzzheng/apps/zMutScan/scripts')

from vespa_scorer import VESPAScorer
from esmif_scorer import ESMIFScorer
from foldx_scorer import FoldXScorer
from base_scorer import Mutation
from utils import pdb2fasta

class UnifiedMutationScorer:
    """Unified interface for multiple scoring methods."""
    
    def __init__(self, pdb_file, output_base="./results"):
        self.pdb_file = pdb_file
        self.output_base = output_base
        self.sequence = pdb2fasta(pdb_file)
        
        # Initialize scorers
        self.vespa_scorer = VESPAScorer(
            config={"args": {"model": "esm2_650M"}},
            output_dir=f"{output_base}/vespa",
            wt_sequence=self.sequence,
            pdb_file=pdb_file
        )
        
        self.esmif_scorer = ESMIFScorer(
            config={"args": {"model": "esm_if1_gvp4_t16_142M_UR50"}},
            output_dir=f"{output_base}/esmif",
            wt_sequence=self.sequence,
            pdb_file=pdb_file
        )
        
        self.foldx_scorer = FoldXScorer(
            config={"args": {"num_runs": 5}},
            output_dir=f"{output_base}/foldx",
            wt_sequence=self.sequence,
            pdb_file=pdb_file
        )
    
    def score_mutation(self, mutation_str):
        """Score mutation with all methods."""
        
        mutation = Mutation.from_string(mutation_str)
        
        scores = {
            "mutation": mutation_str,
            "vespa": self.vespa_scorer.score_mutation(mutation),
            "esmif": self.esmif_scorer.score_mutation(mutation),
            "foldx": self.foldx_scorer.score_mutation(mutation),
        }
        
        # Calculate consensus
        scores["consensus"] = self._calculate_consensus(scores)
        
        return scores
    
    def _calculate_consensus(self, scores):
        """Calculate consensus score from multiple methods."""
        # Normalize scores (all lower = better)
        normalized = [
            scores["vespa"],
            scores["esmif"],
            scores["foldx"]
        ]
        
        # Remove None values
        normalized = [s for s in normalized if s is not None]
        
        if normalized:
            return sum(normalized) / len(normalized)
        return None
    
    def batch_score(self, mutations):
        """Score multiple mutations."""
        results = []
        for mut in mutations:
            try:
                scores = self.score_mutation(mut)
                results.append(scores)
                print(f"Completed: {mut}")
            except Exception as e:
                print(f"Error for {mut}: {e}")
                results.append({"mutation": mut, "error": str(e)})
        
        return results

# Example usage
if __name__ == "__main__":
    scorer = UnifiedMutationScorer(
        pdb_file="protein.pdb",
        output_base="./unified_results"
    )
    
    mutations = ["M198F", "K201E", "D205N"]
    results = scorer.batch_score(mutations)
    
    for result in results:
        print(f"\n{result['mutation']}:")
        print(f"  VESPA: {result.get('vespa', 'N/A'):.6f}")
        print(f"  ESM-IF: {result.get('esmif', 'N/A'):.6f}")
        print(f"  FoldX: {result.get('foldx', 'N/A'):.3f}")
        print(f"  Consensus: {result.get('consensus', 'N/A'):.6f}")
```

## Score Interpretation Guide

### VESPA Scores
- **< -2**: Strongly favorable
- **-2 to 0**: Favorable to neutral
- **0 to 2**: Unfavorable
- **> 2**: Strongly unfavorable

### ESM-IF Scores (ΔLoss)
- **< -0.5**: Strongly favorable (better structure fit)
- **-0.5 to -0.1**: Favorable
- **-0.1 to 0.1**: Neutral
- **0.1 to 0.5**: Unfavorable
- **> 0.5**: Strongly unfavorable (worse structure fit)

### FoldX Scores (ΔΔG in kcal/mol)
- **< -1.0**: Strongly stabilizing
- **-1.0 to -0.5**: Stabilizing
- **-0.5 to 0.5**: Neutral
- **0.5 to 1.0**: Destabilizing
- **> 1.0**: Strongly destabilizing

### Consensus Interpretation
```python
def interpret_consensus(scores):
    """Interpret consensus from multiple scoring methods."""
    
    favorable_count = 0
    total_count = 0
    
    # VESPA: < 0 is favorable
    if scores.get("vespa") is not None:
        total_count += 1
        if scores["vespa"] < 0:
            favorable_count += 1
    
    # ESM-IF: < 0 is favorable
    if scores.get("esmif") is not None:
        total_count += 1
        if scores["esmif"] < 0:
            favorable_count += 1
    
    # FoldX: < 0 is favorable
    if scores.get("foldx") is not None:
        total_count += 1
        if scores["foldx"] < 0:
            favorable_count += 1
    
    if total_count == 0:
        return "No scores available"
    
    favorable_ratio = favorable_count / total_count
    
    if favorable_ratio >= 0.67:
        return "Consensus favorable"
    elif favorable_ratio >= 0.33:
        return "Mixed evidence"
    else:
        return "Consensus unfavorable"
```

## Performance Optimization

### Parallel Processing
```python
#!/usr/bin/env python3
"""
Parallel mutation scoring
"""

from concurrent.futures import ProcessPoolExecutor, as_completed

def parallel_score_mutations(pdb_file, mutations, output_base, max_workers=4):
    """Score mutations in parallel."""
    
    def score_single(mut):
        """Score single mutation in separate process."""
        import sys
        sys.path.insert(0, '/data_test/home/lzzheng/apps/zMutScan/scripts')
        
        from unified_scorer import UnifiedMutationScorer
        scorer = UnifiedMutationScorer(pdb_file, f"{output_base}/{mut}")
        return scorer.score_mutation(mut)
    
    results = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_mut = {executor.submit(score_single, mut): mut for mut in mutations}
        
        for future in as_completed(future_to_mut):
            mut = future_to_mut[future]
            try:
                results[mut] = future.result()
                print(f"✓ Completed: {mut}")
            except Exception as e:
                results[mut] = {"error": str(e)}
                print(f"✗ Failed: {mut} - {e}")
    
    return results
```

### Caching Strategy
```python
import json
import os
from datetime import datetime

class ScoreCache:
    """Cache for mutation scores."""
    
    def __init__(self, cache_file=".mutation_scores_cache.json"):
        self.cache_file = cache_file
        self.cache = self._load_cache()
    
    def _load_cache(self):
        """Load cache from file."""
        if os.path.exists(self.cache_file):
            with open(self.cache_file, 'r') as f:
                return json.load(f)
        return {}
    
    def save_cache(self):
        """Save cache to file."""
        with open(self.cache_file, 'w') as f:
            json.dump(self.cache, f, indent=2)
    
    def get(self, pdb_hash, mutation, tool):
        """Get cached score."""
        key = f"{pdb_hash}:{mutation}:{tool}"
        return self.cache.get(key)
    
    def set(self, pdb_hash, mutation, tool, score):
        """Set cached score."""
        key = f"{pdb_hash}:{mutation}:{tool}"
        self.cache[key] = {
            "score": score,
            "timestamp": datetime.now().isoformat(),
            "pdb": pdb_hash,
            "mutation": mutation,
            "tool": tool
        }
        self.save_cache()
    
    def clear_old_entries(self, days=30):
        """Clear entries older than specified days."""
        cutoff = datetime.now().timestamp() - (days * 24 * 3600)
        
        new_cache = {}
        for key, entry in self.cache.items():
            entry_time = datetime.fromisoformat(entry["timestamp"]).timestamp()
            if entry_time > cutoff:
                new_cache[key] = entry
        
        self.cache = new_cache
        self.save_cache()
```

## Visualization

### Multi-Tool Score Plot
```python
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_mutation_scores(results_df, output_file="mutation_scores.png"):
    """Plot mutation scores from multiple tools."""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # VESPA scores
    ax1 = axes[0, 0]
    vespa_data = results_df[["mutation", "vespa"]].dropna()
    ax1.bar(range(len(vespa_data)), vespa_data["vespa"])
    ax1.set_xticks(range(len(vespa_data)))
    ax1.set_xticklabels(vespa_data["mutation"], rotation=45, ha='right')
    ax1.set_title("VESPA Scores")
    ax1.set_ylabel("Score (lower = better)")
    ax1.axhline(y=0, color='r', linestyle='-', alpha=0.3)
    
    # ESM-IF scores
    ax2 = axes[0, 1]
    esmif_data = results_df[["mutation", "esmif"]].dropna()
    ax2.bar(range(len(esmif_data)), esmif_data["esmif"])
    ax2.set_xticks(range(len(esmif_data)))
    ax2.set_xticklabels(esmif_data["mutation"], rotation=45, ha='right')
    ax2.set_title("ESM-IF Scores")
    ax2.set_ylabel("ΔLoss (lower = better)")
    ax2.axhline(y=0, color='r', linestyle='-', alpha=0.3)
    
    # FoldX scores
    ax3 = axes[1, 0]
    foldx_data = results_df[["mutation", "foldx"]].dropna()
    ax3.bar(range(len(foldx_data)), foldx_data["foldx"])
    ax3.set_xticks(range(len(foldx_data)))
    ax3.set_xticklabels(foldx_data["mutation"], rotation=45, ha='right')
    ax3.set_title("FoldX Scores")
    ax3.set_ylabel("ΔΔG (kcal/mol, lower = better)")
    ax3.axhline(y=0, color='r', linestyle='-', alpha=0.3)
    
    # Consensus scores
    ax4 = axes[1, 1]
    consensus_data = results_df[["mutation", "consensus"]].dropna()
    ax4.bar(range(len(consensus_data)), consensus_data["consensus"])
    ax4.set_xticks(range(len(consensus_data)))
    ax4.set_xticklabels(consensus_data["mutation"], rotation=45, ha='right')
    ax4.set_title("Consensus Scores")
    ax4.set_ylabel("Normalized Score")
    ax4.axhline(y=0, color='r', linestyle='-', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
```

## Troubleshooting Common Issues

### 1. Tool Not Found
```bash
# Check tool availability
which FoldX  # FoldX binary
python -c "import esm; print(esm.__version__)"  # ESM package
```

### 2. Memory Issues
```bash
# Reduce batch size
export VESPA_BATCH_SIZE=1
export FOLDX_NUM_RUNS=3

# Use smaller models
python script.py --model esm2_150M  # Instead of esm2_3B
```

### 3. PDB Format Problems
```python
# Clean PDB file
from Bio.PDB import PDBParser, PDBIO
parser = PDBParser()
structure = parser.get_structure("clean", "protein.pdb")
io = PDBIO()
io.set_structure(structure)
io.save("protein_clean.pdb")
```

### 4. Score Interpretation
```python
# Always check score ranges
print(f"VESPA range: {df['vespa'].min():.3f} to {df['vespa'].max():.3f}")
print(f"ESM-IF range: {df['esmif'].min():.3f} to {df['esmif'].max():.3f}")
print(f"FoldX range: {df['foldx'].min():.3f} to {df['foldx'].max():.3f}")
```

## Best Practices

1. **Use multiple tools**: No single tool is perfect
2. **Validate with experimental data**: When available
3. **Clean PDB files**: Essential for structure-based methods
4. **Consider biological context**: Active sites, interfaces, etc.
5. **Document parameters**: Record model versions and settings

## References

- **VESPA**: Marquet et al., Human Genetics 2021
- **ESM-IF**: Hsu et al., ICML 2022
- **FoldX**: Schymkowitz et al., NAR 2005
- **ESM2**: Lin et al., bioRxiv 2022
- **zMutScan**: Local mutation scanning framework