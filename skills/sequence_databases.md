# Protein Sequence Databases Skill

## Quick Start

### Generate PSSM from FASTA
```bash
# 1. Generate MSA
qjackhmmer -N 3 -B query.a3m --incE 0.001 query.fasta /data_test/share/pub_datas/af3/uniref90_2022_05.fa > /dev/null 2>&1

# 2. Convert to PSSM
MSA_To_PSSM -I query.a3m -o query.pssm -C 32

# 3. Parse PSSM in Python
from zSaprod.scorer.pssm import PSSMParser
parser = PSSMParser("query.pssm")
```

## Overview
This skill provides access to protein sequence databases for Multiple Sequence Alignment (MSA) generation and Position-Specific Scoring Matrix (PSSM) calculation. The databases are located at `/data_test/share/pub_datas/af3/` and can be used with tools like `qjackhmmer` for MSA generation and `MSA_To_PSSM` for PSSM calculation.

## Available Databases

### Protein Databases
1. **uniprot_all_2021_04.fa** (84GB) - Complete UniProt database (April 2021)
   - Comprehensive protein sequences from UniProt
   - Best for general protein homology searches

2. **uniref90_2022_05.fa** (67GB) - UniRef90 clustered database (May 2022)
   - Clustered at 90% sequence identity
   - Reduces redundancy while maintaining diversity
   - Faster searches than full UniProt

3. **pdb_seqres_2022_09_28.fasta** (223MB) - PDB sequences (September 2022)
   - Sequences from Protein Data Bank structures
   - Useful for structure-based homology

4. **mgy_clusters_2022_05.fa** (120GB) - Metagenomic clusters (May 2022)
   - Metagenomic protein sequences
   - Useful for novel or environmental proteins

5. **bfd-first_non_consensus_sequences.fasta** (17GB) - BFD non-consensus sequences
   - Big Fantastic Database (BFD) sequences
   - Used in AlphaFold and other structure prediction tools

### RNA/DNA Databases
1. **nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta** (76GB) - Nucleotide/RNA database
   - Clustered at 90% identity, 80% coverage
   - For RNA/protein-RNA interaction studies

2. **rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta** (218MB) - Rfam RNA families
   - RNA family database
   - For non-coding RNA homology

3. **rnacentral_active_seq_id_90_cov_80_linclust.fasta** (13GB) - RNAcentral active sequences
   - Comprehensive non-coding RNA sequences

## Tools

### qjackhmmer (Iterative HMMER Search)
Location: `~/bin/qjackhmmer` or `/data_test/home/lzzheng/apps/autoFuncLib/tools/qjackhmmer`

```bash
# Basic usage
qjackhmmer -N 3 -B output.a3m --incE 0.001 query.fasta /data_test/share/pub_datas/af3/uniprot_all_2021_04.fa

# Common options:
# -N <n>    : maximum iterations (default: 5)
# -B <file> : save alignment to A3M file
# --incE <x>: inclusion E-value threshold
# --cpu <n> : number of CPU cores to use
```

### MSA_To_PSSM (Convert MSA to PSSM)
Location: `~/bin/MSA_To_PSSM` or `/data_test/home/lzzheng/apps/autoFuncLib/tools/MSA_To_PSSM`

### MSA_CovFilter (Filter MSA by coverage)
Location: `/data_test/home/lzzheng/apps/autoFuncLib/tools/MSA_CovFilter`

```bash
# Filter MSA by coverage threshold
MSA_CovFilter input.a3m output.a3m 80

# 80 = 80% coverage threshold
```

```bash
# Basic usage
MSA_To_PSSM -I input.a3m -o output.pssm -C 32 -c 10000

# Options:
# -I <file> : input A3M/PSI format MSA file
# -o <file> : output PSSM file
# -C <n>    : CPU number (default: all)
# -c <n>    : CD-HIT cutoff for large MSAs
```

## Database Sources

These databases were fetched from Google DeepMind's AlphaFold 3 database repository using the `fetch_databases.sh` script. They are the same databases used by AlphaFold 3 for structure prediction.

### Database Fetch Script
```bash
# Script location: /data_test/share/pub_datas/af3/fetch_databases.sh
# Source: https://storage.googleapis.com/alphafold-databases/v3.0
```

## Python PSSM Parser

The PSSM parser is available at: `/data_test/home/lzzheng/apps/zSaprod/zSaprod/scorer/pssm.py`

### Basic Usage
```python
from zSaprod.scorer.pssm import PSSMParser, PSSMScorer, GeneratePSSM

# Parse PSSM file
parser = PSSMParser("result.pssm")
print("Sequence:", "".join(parser._aa_sequence))
print("PSSM matrix shape:", len(parser._lpssmc), "x", len(parser._lpssmc[0]))

# Score a sequence
scorer = PSSMScorer(pssm_fpath="result.pssm")
score = scorer.score(sequence_object)
```

### Generate PSSM from FASTA
```python
from zSaprod.scorer.pssm import GeneratePSSM
import json

# Load configuration
with open("config.json") as f:
    configs = json.load(f)

# Generate PSSM
gpssm = GeneratePSSM(
    fasta_fpath="query.fasta",
    output_dpath="output/pssm",
    configs=configs['scorer']['pssm']
)
gpssm.pssm_generation()
```

## Workflow Examples

### 1. Generate MSA using qjackhmmer
```bash
# Using UniRef90 (faster, less redundant)
qjackhmmer -N 3 -B query.a3m --incE 0.001 query.fasta /data_test/share/pub_datas/af3/uniref90_2022_05.fa > /dev/null 2>&1

# Using full UniProt (more comprehensive)
qjackhmmer -N 3 -B query.a3m --incE 0.001 query.fasta /data_test/share/pub_datas/af3/uniprot_all_2021_04.fa > /dev/null 2>&1

# Using multiple CPUs
qjackhmmer -N 3 -B query.a3m --incE 0.001 --cpu 32 query.fasta /data_test/share/pub_datas/af3/uniprot_all_2021_04.fa > /dev/null 2>&1
```

### 2. Convert MSA to PSSM
```bash
# Convert A3M to PSSM
MSA_To_PSSM -I query.a3m -o query.pssm -C 32

# With CD-HIT filtering for large MSAs
MSA_To_PSSM -I query.a3m -o query.pssm -C 32 -c 10000
```

### 3. Complete Pipeline Script
```bash
#!/bin/bash
# complete_pssm_pipeline.sh

QUERY_FASTA=$1
OUTPUT_DIR=$2
DATABASE=$3

# Create output directory
mkdir -p $OUTPUT_DIR

# Step 1: Generate MSA
echo "Generating MSA..."
qjackhmmer -N 3 -B $OUTPUT_DIR/query.a3m --incE 0.001 $QUERY_FASTA $DATABASE > /dev/null 2>&1

# Step 2: Convert to PSSM
echo "Converting MSA to PSSM..."
MSA_To_PSSM -I $OUTPUT_DIR/query.a3m -o $OUTPUT_DIR/query.pssm -C 32

echo "PSSM generated at: $OUTPUT_DIR/query.pssm"
```

### 4. Python Integration
```python
import subprocess
import os

def generate_pssm_pipeline(query_fasta, output_dir, database_path):
    """Complete PSSM generation pipeline."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate MSA
    a3m_file = os.path.join(output_dir, "query.a3m")
    cmd = f"qjackhmmer -N 3 -B {a3m_file} --incE 0.001 {query_fasta} {database_path}"
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # Convert to PSSM
    pssm_file = os.path.join(output_dir, "query.pssm")
    cmd = f"MSA_To_PSSM -I {a3m_file} -o {pssm_file} -C 32"
    subprocess.run(cmd, shell=True)
    
    return pssm_file
```

## Database Selection Guide

| Use Case | Recommended Database | Reason |
|----------|---------------------|--------|
| General protein homology | `uniref90_2022_05.fa` | Fast, less redundant, good coverage |
| Comprehensive search | `uniprot_all_2021_04.fa` | Most comprehensive, includes all UniProt |
| Structure-based homology | `pdb_seqres_2022_09_28.fasta` | Only sequences with known structures |
| Novel/environmental proteins | `mgy_clusters_2022_05.fa` | Metagenomic sequences |
| AlphaFold-style MSA | `bfd-first_non_consensus_sequences.fasta` | Used in AlphaFold training |
| RNA/protein-RNA | `nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta` | RNA sequences |

## Performance Tips

1. **Use UniRef90 for speed**: 67GB vs 84GB for full UniProt
2. **Limit iterations**: `-N 3` is usually sufficient
3. **Use multiple CPUs**: `--cpu 32` for parallel processing
4. **Filter MSA**: Use `-c 10000` with `MSA_To_PSSM` for large alignments
5. **Redirect output**: Use `> /dev/null 2>&1` to suppress qjackhmmer output

## Configuration File Example

Create a `config.json` file for the PSSM generator (based on actual configs_wzh.json):

```json
{
  "scorer": {
    "pssm": {
      "qjackhmm": "/data_test/home/lzzheng/apps/autoFuncLib/tools/qjackhmmer",
      "MSA2PSSM": "/data_test/home/lzzheng/apps/autoFuncLib/tools/MSA_To_PSSM",
      "msa_filter": "/data_test/home/lzzheng/apps/autoFuncLib/tools/MSA_CovFilter",
      "ncpus": 16,
      "db": "/data_test/share/pub_datas/af3/uniref90_2022_05.fa",
      "result_fpath": "result.pssm"
    }
  }
}
```

## Troubleshooting

### Common Issues

1. **qjackhmmer not found**: Ensure `~/bin/` is in your PATH
2. **Database path errors**: Use absolute paths to databases
3. **Memory issues**: For large queries, use smaller database (UniRef90 instead of full UniProt)
4. **Slow performance**: Use `--cpu` option with available cores

### Checking Database Integrity
```bash
# Check if database is accessible
head -5 /data_test/share/pub_datas/af3/uniref90_2022_05.fa

# Check file size
ls -lh /data_test/share/pub_datas/af3/uniref90_2022_05.fa
```

## References

- **HMMER**: http://hmmer.org/
- **UniProt**: https://www.uniprot.org/
- **UniRef**: https://www.uniprot.org/help/uniref
- **PDB**: https://www.rcsb.org/
- **BFD**: https://bfd.mmseqs.com/