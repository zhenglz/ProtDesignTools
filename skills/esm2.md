# ESM2 Skill

## Overview
ESM2 (Evolutionary Scale Modeling 2) is a transformer-based protein language model that learns evolutionary patterns from protein sequences. It can be used for sequence embeddings, variant effect prediction, and structure prediction.

## Installation
```bash
# Install via pip
pip install fair-esm

# Or install from source
git clone https://github.com/facebookresearch/esm.git
cd esm
pip install -e .
```

## Basic Usage

### Command Line Interface
```bash
# Generate embeddings for a FASTA file
python scripts/extract.py esm2_t33_650M_UR50D \
  <input_fasta> \
  <output_directory> \
  --repr_layers 33 \
  --include mean

# Example
python scripts/extract.py esm2_t33_650M_UR50D \
  sequences.fasta \
  ./embeddings \
  --repr_layers 33 \
  --include mean per_tok
```

### Key Parameters
- `model_name`: ESM2 model variant (esm2_t6_8M_UR50D to esm2_t48_15B_UR50D)
- `input_fasta`: Input FASTA file
- `output_dir`: Output directory for embeddings
- `--repr_layers`: Which layers to extract (comma-separated or "all")
- `--include`: What to include ("mean", "per_tok", "bos", "contacts")
- `--truncation_seq_length`: Maximum sequence length (default: 1024)
- `--nogpu`: Disable GPU (use CPU)

## Input Requirements

### Input Files
1. **FASTA file**: Standard protein sequence file
2. **Single sequence**: Can also pass sequence directly via stdin

### Input Preparation
```bash
# Ensure FASTA format
cat > sequences.fasta << 'EOF'
>protein1
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLA
>protein2
MKLKKLRLPSLKVLLFIFSLLLACSSPSAQTSPAVQTSPAVQTSPAVQTSPAVQTSPAVQTSP
EOF

# Check sequence length
python -c "
from Bio import SeqIO
for record in SeqIO.parse('sequences.fasta', 'fasta'):
    print(f'{record.id}: {len(record.seq)} residues')
"
```

## Output Structure

### Generated Files
```
output_directory/
├── <sequence_id>_layer_<layer_number>.pt    # Per-token embeddings (PyTorch tensors)
├── <sequence_id>_mean_representations.pt    # Mean embeddings
├── <sequence_id>_contacts.pt                # Contact maps
└── <sequence_id>_log.txt                    # Processing log
```

### Output Files Description

1. **PyTorch tensor files** (`*.pt`):
   - Per-token embeddings: `[sequence_length, embedding_dim]`
   - Mean embeddings: `[embedding_dim]`
   - Contact maps: `[sequence_length, sequence_length]`

2. **Text format** (optional):
   - Can output as numpy `.npy` or CSV format
   - Use `--toks_per_batch` for batching long sequences

## Output Interpretation

### Embedding Analysis
```bash
# Load and analyze embeddings
python -c "
import torch
import numpy as np

# Load embeddings
embeddings = torch.load('embeddings/protein1_layer_33.pt')
print(f'Embedding shape: {embeddings.shape}')
print(f'Mean embedding norm: {torch.norm(embeddings, dim=1).mean():.2f}')

# Calculate sequence similarity
if len(embeddings.shape) == 2:  # mean embeddings
    emb1 = torch.load('embeddings/protein1_mean_representations.pt')
    emb2 = torch.load('embeddings/protein2_mean_representations.pt')
    similarity = torch.cosine_similarity(emb1, emb2, dim=0)
    print(f'Cosine similarity: {similarity:.3f}')
"
```

### Variant Effect Prediction
```bash
# Predict variant effects
python scripts/variant_prediction.py \
  --model-location esm2_t33_650M_UR50D \
  --sequence MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLA \
  --dms-input variants.csv \
  --mutation-col mutant \
  --dms-output variant_scores.csv \
  --offset-idx 1 \
  --scoring-strategy wt-marginals
```

## Advanced Usage

### Multiple Layers and Representations
```bash
# Extract multiple layers and representations
python scripts/extract.py esm2_t33_650M_UR50D \
  sequences.fasta \
  ./embeddings \
  --repr_layers 0,6,12,18,24,30,33 \
  --include mean per_tok contacts \
  --toks_per_batch 1000
```

### Structure Prediction (ESMFold)
```bash
# Predict structure with ESMFold
python scripts/esmfold.py \
  --model esmfold_v1 \
  --fasta sequences.fasta \
  --output-dir ./structures \
  --num-recycles 4 \
  --chunk-size 128

# Or use the API
import esm
model = esm.pretrained.esmfold_v1()
model = model.eval().cuda()
output = model.infer("MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLA")
pdbs = output["pdbs"]
```

### Contact Map Prediction
```bash
# Predict contact maps
python scripts/extract.py esm2_t33_650M_UR50D \
  sequences.fasta \
  ./contacts \
  --include contacts \
  --filter-threshold 0.2

# Visualize contacts
python scripts/visualize_contacts.py \
  --contacts contacts/protein1_contacts.pt \
  --sequence sequences.fasta \
  --output contact_map.png
```

## Troubleshooting

### Common Issues
1. **CUDA out of memory**: Reduce batch size with `--toks_per_batch`
2. **Sequence too long**: Use `--truncation_seq_length` or split sequences
3. **Model download errors**: Check internet connection and disk space

### Debug Commands
```bash
# Test with small sequence
echo ">test\nMKTVRQERLK" > test.fasta
python scripts/extract.py esm2_t6_8M_UR50D test.fasta ./test_out

# Check GPU availability
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"

# List available models
python -c "import esm; print(esm.pretrained.available_models())"
```

## Integration with Other Tools

### Downstream Analysis
```bash
# Cluster sequences by embeddings
python scripts/cluster_embeddings.py \
  --embeddings-dir ./embeddings \
  --output clusters.csv \
  --method kmeans \
  --n-clusters 10

# Train classifier on embeddings
python scripts/train_classifier.py \
  --embeddings ./embeddings \
  --labels labels.csv \
  --output-classifier classifier.pkl
```

## Model Variants

| Model | Parameters | Layers | Embedding Dim | Context |
|-------|------------|--------|---------------|---------|
| esm2_t6_8M_UR50D | 8M | 6 | 320 | 1024 |
| esm2_t12_35M_UR50D | 35M | 12 | 480 | 1024 |
| esm2_t30_150M_UR50D | 150M | 30 | 640 | 1024 |
| esm2_t33_650M_UR50D | 650M | 33 | 1280 | 1024 |
| esm2_t36_3B_UR50D | 3B | 36 | 2560 | 1024 |
| esm2_t48_15B_UR50D | 15B | 48 | 5120 | 1024 |

## References
- [ESM GitHub](https://github.com/facebookresearch/esm)
- [ESM2 Paper](https://www.science.org/doi/10.1126/science.ade2574)
- [ESMFold Paper](https://www.science.org/doi/10.1126/science.ade2574)