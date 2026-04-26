# Literature Miner Tool

**File:** `tools/literature_miner.py`

A comprehensive tool for searching, downloading, and querying academic literature. Supports three modes: web search & PDF download, local PDF ingestion, and RAG-based Q&A over a cached knowledge database.

---

## Table of Contents

- [Overview](#overview)
- [Requirements & Setup](#requirements--setup)
- [Modes of Operation](#modes-of-operation)
  - [Mode 1: Web Search & Download (`--query`)](#mode-1-web-search--download---query)
  - [Mode 2: PDF Directory Ingestion (`--pdf-dir`)](#mode-2-pdf-directory-ingestion---pdf-dir)
  - [Mode 3: Database Search (`--query-db`)](#mode-3-database-search---query-db)
  - [Mode 4: Agent RAG Q&A (`--agent`)](#mode-4-agent-rag-qa---agent)
- [Knowledge Database](#knowledge-database)
- [Output Formats](#output-formats)
- [Architecture & Data Flow](#architecture--data-flow)
  - [Web Search Pipeline](#web-search-pipeline)
  - [PDF Ingestion Pipeline](#pdf-ingestion-pipeline)
  - [Database Caching Strategy](#database-caching-strategy)
  - [Search & Retrieval](#search--retrieval)
- [Configuration](#configuration)
- [Full Argument Reference](#full-argument-reference)
- [Examples](#examples)

---

## Overview

The Literature Miner tool automates the process of:

1. **Discovering** academic papers via web search (Tavily API)
2. **Downloading** PDFs from arXiv, PubMed Central, Europe PMC, Sci-Hub, and direct publisher links
3. **Extracting** text from PDFs using PyPDF2 / pypdf / pdfplumber
4. **Summarizing** papers using an LLM (DeepSeek by default, via OpenAI-compatible API)
5. **Storing** structured metadata in a local JSON knowledge database
6. **Querying** the database with keyword search or natural-language RAG Q&A

The knowledge database is JSON-based and acts as a persistent cache — re-running the tool only processes new/changed PDFs.

---

## Requirements & Setup

### Python Dependencies

| Dependency | Required For | Install Command |
|---|---|---|
| `requests` | API calls (Tavily, DeepSeek) | `pip install requests` |
| `pypdf` | PDF text extraction | `pip install pypdf` |
| `PyPDF2` | PDF extraction (alternative) | `pip install PyPDF2` |
| `pdfplumber` | PDF extraction (alternative) | `pip install pdfplumber` |

At least one PDF extraction library is required for `--pdf-dir` mode and for the summarization step in `--query` mode.

### API Keys

| Key | Environment Variable | CLI Flag | Source |
|---|---|---|---|
| Tavily | `TAVILY_API_KEY` | `--tavily-api-key` | [tavily.com](https://tavily.com) |
| DeepSeek | `OPENAI_API_KEY` | `--deepseek-api-key` | [platform.deepseek.com](https://platform.deepseek.com) |

The tool automatically loads a `.env` file from the repository root (searched upward from the script's location). The `.env` file uses bash `export` format:

```bash
export TAVILY_API_KEY="tvly-..."
export OPENAI_API_KEY="sk-..."
export OPENAI_BASE_URL="https://api.deepseek.com/v1"
export OPENAI_MODEL="deepseek-chat"
```

Keys provided via CLI flags take precedence over environment variables.

---

## Modes of Operation

The four modes are **mutually exclusive** — specify exactly one.

### Mode 1: Web Search & Download (`--query`)

Search the web for academic papers, identify them via DOIs / arXiv IDs / PubMed IDs, download PDFs, summarize with LLM, and add to the knowledge database.

```
python tools/literature_miner.py --query "protein diffusion models" \
  --max-downloads 5 --output-dir ./literature_mining
```

**Pipeline:**
1. **Phase 1 — Web Search**: Issues the query to Tavily, plus site-specific queries for open-access sources (`arxiv.org`, `biorxiv.org`, `frontiersin.org`, `mdpi.com`, `plos.org`). Deduplicates by URL.
2. **Phase 2 — Paper Identification**: Extracts DOIs, arXiv IDs, PubMed IDs, PMC IDs from each result URL and content snippet. Filters to results that have at least one academic identifier or come from an academic domain.
3. **Phase 3 — PDF Download**: Builds prioritized download URL lists for each paper (arXiv > Europe PMC > PubMed Central > Sci-Hub mirrors > page scraping). Validates responses by checking PDF magic bytes (`%PDF-`).
4. **Phase 4 — LLM Summarization**: Extracts PDF text (up to 8000 chars), sends to DeepSeek with a structured JSON prompt requesting summary, key findings, and keywords.
5. **Phase 5 — Database Storage**: Deduplicates by DOI/arXiv ID/PubMed ID, records search history.
6. **Phase 6 — Report**: Prints results in the chosen output format.

**Flags:**
- `--no-download`: Skip PDF download, collect metadata only
- `--no-summarize`: Skip LLM summarization (no DeepSeek key needed)
- `--max-results N`: Max search results to process (default: 20)
- `--max-downloads N`: Max PDFs to download (default: 10)
- `--scihub-mirror URL`: Sci-Hub mirror URL (default: `https://sci-hub.se`)

### Mode 2: PDF Directory Ingestion (`--pdf-dir`)

Read PDF files from a local directory, extract text, summarize with DeepSeek, and add to the knowledge database. **Caches results by SHA-256 file hash** — re-running skips unchanged PDFs.

```
python tools/literature_miner.py --pdf-dir ./my_papers \
  --db-path ./literature_mining/knowledge_db.json
```

**Pipeline:**
1. Loads the existing knowledge database (or creates a new one)
2. Builds a set of known SHA-256 file hashes from existing entries
3. Scans the directory for `*.pdf` files (case-insensitive, sorted alphabetically)
4. For each PDF:
   - Computes `hash_file()` — SHA-256 hash of the full file content
   - Skips if the hash exists in the database and `--pdf-skip-existing` is enabled (default)
   - Extracts text via `extract_pdf_text()` (tries PyPDF2 → pdfplumber → pypdf)
   - Sends text to DeepSeek for summarization
   - Stores a structured entry with `file_hash`, summary, key findings, keywords
5. Saves the updated database

**Flags:**
- `--pdf-max-chars N`: Max characters to extract from each PDF (default: 8000)
- `--no-pdf-skip-existing`: Re-process all PDFs regardless of cache status

### Mode 3: Database Search (`--query-db`)

Search the local knowledge database for papers matching a keyword or phrase.

```
python tools/literature_miner.py --query-db "diffusion" \
  --db-path ./literature_mining/knowledge_db.json
```

**Search algorithm:**
- Splits the query into terms
- For each paper, scores match across fields with different weights:
  - Title: weight 10, boost ×5 for exact phrase match
  - Keywords: weight 8, boost ×5
  - Summary: weight 6, boost ×5
  - Key findings: weight 5, boost ×5
  - Abstract: weight 3, boost ×5
  - Original search query: weight 2, boost ×5
- Returns papers sorted by descending score

### Mode 4: Agent RAG Q&A (`--agent`)

Ask natural-language questions against the knowledge database. Uses DeepSeek to retrieve relevant papers and generate a synthesized answer with citations.

**Single question:**
```
python tools/literature_miner.py --agent \
  --question "What are key advances in protein diffusion?" \
  --db-path ./literature_mining/knowledge_db.json
```

**Interactive session:**
```
python tools/literature_miner.py --agent --agent-interactive \
  --db-path ./literature_mining/knowledge_db.json
```

**Interactive commands:**
| Command | Action |
|---|---|
| `exit` / `quit` | End session |
| `help` | Show available commands |
| `stats` | Show paper count and unique keywords |
| `papers` | List all paper titles |
| *any question* | Run RAG query on current database |

**How it works:**
1. Retrieves top-K relevant papers using `search_db()` (weighted keyword scoring)
2. Builds a context string with title, authors, summary, findings, and links for each paper
3. Sends context + question to DeepSeek with a system prompt instructing it to cite papers using `[Paper N]` notation
4. Returns the answer along with source papers

**Flags:**
- `--question TEXT`: The question to answer
- `--agent-top-k N`: Number of papers to retrieve for context (default: 5)
- `--agent-interactive`: Run in interactive REPL mode

---

## Knowledge Database

The database is a JSON file with this structure:

```json
{
  "papers": [
    {
      "id": "a1b2c3d4e5f6",
      "title": "Protein Diffusion Models: A Comprehensive Review",
      "authors": ["John Smith", "Jane Doe"],
      "abstract": "Diffusion models have emerged as...",
      "source_url": "https://arxiv.org/abs/1234.56789",
      "doi": "10.1234/example.doi",
      "arxiv_id": "1234.56789",
      "pubmed_id": null,
      "pmcid": null,
      "pdf_path": "/path/to/paper.pdf",
      "file_hash": "sha256hex...",
      "keywords": ["diffusion models", "protein design", "generative models"],
      "summary": "This review surveys recent advances...",
      "key_findings": [
        "Diffusion models are a powerful class of generative models...",
        "Key methods include DDPM, score-based generative models..."
      ],
      "search_query": "protein diffusion models",
      "date_added": "2026-04-26"
    }
  ],
  "search_history": [
    {
      "query": "protein diffusion models",
      "date": "2026-04-26 22:30:00",
      "results_count": 15,
      "new_papers": 3
    }
  ]
}
```

**Deduplication:** Papers are kept unique by DOI, arXiv ID, and PubMed ID. The `add_paper_to_db()` function checks these identifiers before inserting.

**For PDF-ingested entries:**
- `doi`, `arxiv_id`, `pubmed_id`, `pmcid` are all `null`
- `source_url` and `pdf_path` point to the local file
- `file_hash` stores the SHA-256 hash for caching
- `title` is the filename (e.g., `protein_diffusion_review.pdf`)
- Authors and abstract are empty (the LLM summary substitutes for them)

---

## Output Formats

Results can be displayed in three formats, controlled by `--output-format`:

| Format | Flag | Description |
|---|---|---|
| Table | `--output-format table` | Formatted text table with paper number, title, DOI, arXiv, source, PDF path, abstract excerpt, summary, findings |
| JSON | `--output-format json` | Raw JSON array of paper objects (excluding score) |
| CSV | `--output-format csv` | CSV with columns: title, doi, arxiv_id, pubmed_id, source_url, pdf_path, abstract, summary, keywords |

---

## Architecture & Data Flow

### Web Search Pipeline

```
User Query
    │
    ▼
┌─────────────────────┐
│  Tavily API Search   │  Multiple queries: bare + site-specific
│  (--max-results)     │  (arxiv.org, biorxiv.org, frontiersin.org, etc.)
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│  Paper ID Extraction │  DOI / arXiv / PubMed / PMCID via regex
│  Academic URL Filter │  Domain whitelist + path pattern match
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│  PDF Download        │  Prioritized URL chain:
│  (--max-downloads)   │  arXiv > Europe PMC > PubMed Central
└────────┬────────────┘  > Sci-Hub (3 mirrors) > Page scraping
         │
         ▼
┌─────────────────────┐
│  Text Extraction     │  PyPDF2 → pdfplumber → pypdf
│  (--max-chars)       │  First N chars of concatenated pages
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│  LLM Summarization   │  DeepSeek: structured JSON output
│  (if --no-summarize) │  {summary, key_findings, keywords}
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│  Knowledge Database  │  Dedup by DOI/arXiv/PubMed → append
│  (knowledge_db.json) │  Record in search_history
└─────────────────────┘
```

### PDF Ingestion Pipeline

```
PDF Directory
    │
    ▼
┌─────────────────────┐
│  Scan for *.pdf      │  Sorted alphabetically
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│  Compute SHA-256     │  hash_file() — 64KB chunked
│  Check Cache         │  file_hash in database → SKIP
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│  Text Extraction     │  pypdf / PyPDF2 / pdfplumber
│  (--pdf-max-chars)   │
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│  LLM Summarization   │  DeepSeek with filename as title
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│  Append to Database  │  Entry with file_hash, summary, etc.
│  (cached for reuse)  │
└─────────────────────┘
```

### Database Caching Strategy

The caching mechanism uses SHA-256 file hashing:

1. On first run, each PDF's `hash_file()` result is stored in the `file_hash` field of its database entry
2. On subsequent runs, existing hashes are loaded into a set
3. Before processing each PDF, its current hash is computed and checked against the set
4. If a match is found AND `--pdf-skip-existing` is enabled (default), the file is skipped
5. `--no-pdf-skip-existing` disables this check and re-processes every PDF

This is robust to:
- **File renames**: hash doesn't change, cache hit
- **File modifications**: hash changes, re-processed
- **Duplicate files**: same hash → one database entry per unique content
- **Partial runs**: only unprocessed files are handled

### Search & Retrieval

**Database search** uses a weighted term-matching algorithm:

- Fields are weighted by importance (title: 10, keywords: 8, summary: 6, findings: 5, abstract: 3, original query: 2)
- Exact phrase match gets a 5× boost on the base weight
- Individual term matches add the base weight once per term per field
- Results are sorted by descending total score

**Agent RAG** builds on this search:
1. Retrieves top-K papers via `search_db()`
2. Constructs a context string with structured paper metadata
3. Sends to DeepSeek with a system prompt that enforces citation format (`[Paper N]`)
4. The LLM synthesizes an answer by reasoning across the retrieved papers

---

## Configuration

All configuration is through CLI arguments or environment variables. There is no config file.

**Environment variables** (loaded automatically from `.env` in repo root):
| Variable | Purpose |
|---|---|
| `TAVILY_API_KEY` | Tavily web search API key |
| `OPENAI_API_KEY` | DeepSeek / OpenAI-compatible API key |
| `OPENAI_BASE_URL` | API base URL (default: `https://api.deepseek.com/v1`) |
| `OPENAI_MODEL` | Model name (default: `deepseek-chat`) |

---

## Full Argument Reference

```
usage: literature_miner.py [-h]
                           [--query QUERY | --query-db QUERY_DB | --agent | --pdf-dir PDF_DIR]
                           [--output-dir OUTPUT_DIR] [--tavily-api-key TAVILY_API_KEY]
                           [--deepseek-api-key DEEPSEEK_API_KEY]
                           [--max-results MAX_RESULTS] [--max-downloads MAX_DOWNLOADS]
                           [--no-download] [--no-summarize]
                           [--scihub-mirror SCIHUB_MIRROR] [--question QUESTION]
                           [--agent-top-k AGENT_TOP_K] [--agent-interactive]
                           [--db-path DB_PATH] [--pdf-max-chars PDF_MAX_CHARS]
                           [--no-pdf-skip-existing]
                           [--output-format {table,json,csv}]

Mode selection (mutually exclusive):
  --query QUERY, -q QUERY    Topic to search the web for
  --query-db QUERY_DB        Search the knowledge database for a keyword
  --agent                    Enable RAG Q&A mode
  --pdf-dir PDF_DIR          Ingest PDFs from a local directory

Web search options:
  --output-dir DIR            Output directory (default: ./literature_mining)
  --max-results N             Max search results to process (default: 20)
  --max-downloads N           Max PDFs to download (default: 10)
  --no-download               Skip PDF download, metadata only
  --no-summarize              Skip LLM summarization
  --scihub-mirror URL         Sci-Hub mirror URL (default: https://sci-hub.se)

Agent options:
  --question TEXT             Question to answer
  --agent-top-k N             Papers for context (default: 5)
  --agent-interactive         Interactive Q&A session

PDF directory options:
  --pdf-max-chars N           Max chars to extract per PDF (default: 8000)
  --no-pdf-skip-existing      Re-process PDFs even if cached

Output options:
  --output-format {table,json,csv}  Output format (default: table)

Database options:
  --db-path PATH              Path to knowledge database JSON file

API key options:
  --tavily-api-key KEY        Tavily API key
  --deepseek-api-key KEY      DeepSeek API key
```

---

## Examples

### 1. Search, download, and summarize papers about protein diffusion

```bash
python tools/literature_miner.py --query "protein diffusion models" \
  --max-downloads 5 \
  --max-results 30 \
  --output-dir ./literature_mining
```

### 2. Search metadata only (no downloads or summarization)

```bash
python tools/literature_miner.py --query "CRISPR gene editing" \
  --no-download \
  --max-results 50 \
  --output-dir ./literature_mining
```

### 3. Ingest a directory of PDFs with caching

```bash
# First run — processes all PDFs
python tools/literature_miner.py --pdf-dir ~/papers \
  --db-path ./literature_mining/knowledge_db.json

# Second run — skips all unchanged PDFs (cached)
python tools/literature_miner.py --pdf-dir ~/papers \
  --db-path ./literature_mining/knowledge_db.json

# Force re-process all PDFs
python tools/literature_miner.py --pdf-dir ~/papers \
  --db-path ./literature_mining/knowledge_db.json \
  --no-pdf-skip-existing
```

### 4. Search the knowledge database

```bash
python tools/literature_miner.py --query-db "diffusion" \
  --db-path ./literature_mining/knowledge_db.json

# Output as JSON
python tools/literature_miner.py --query-db "CRISPR" \
  --db-path ./literature_mining/knowledge_db.json \
  --output-format json
```

### 5. Single RAG question

```bash
python tools/literature_miner.py --agent \
  --question "Compare diffusion models and AlphaFold for protein structure" \
  --db-path ./literature_mining/knowledge_db.json \
  --agent-top-k 10
```

### 6. Interactive RAG session

```bash
python tools/literature_miner.py --agent --agent-interactive \
  --db-path ./literature_mining/knowledge_db.json
```

Inside the interactive session:
```
Ask> What are the latest advances in enzyme engineering?
Ask> stats
Ask> papers
Ask> How do these compare to traditional directed evolution?
Ask> exit
```

### 7. Ingest PDFs and immediately query

```bash
# Step 1: Ingest
python tools/literature_miner.py --pdf-dir ~/papers \
  --db-path ./literature_mining/knowledge_db.json

# Step 2: Search
python tools/literature_miner.py --query-db "enzyme" \
  --db-path ./literature_mining/knowledge_db.json

# Step 3: Ask
python tools/literature_miner.py --agent \
  --question "Summarize the key findings about enzyme engineering" \
  --db-path ./literature_mining/knowledge_db.json
```
