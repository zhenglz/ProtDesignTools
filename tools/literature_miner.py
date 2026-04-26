#!/usr/bin/env python3
"""
Literature Miner Tool
=====================

Search the web for academic literature given a topic, download papers (PDFs),
and build a local searchable knowledge database using an LLM for summarization.
Includes an agent mode for RAG-based Q&A over the knowledge database.

Three modes of operation:
  - Web Search & Download: --query "topic"
  - Database Search:       --query-db "keyword"
  - Agent RAG Q&A:         --agent --question "your question"

Usage:
  # Search and download papers
  python tools/literature_miner.py --query "protein diffusion models" \\
    --max-downloads 5 --output-dir ./literature_mining

  # Search and skip downloads (metadata only)
  python tools/literature_miner.py --query "CRISPR gene editing" \\
    --no-download

  # Ingest PDFs from a local directory (with caching)
  python tools/literature_miner.py --pdf-dir ./my_papers \\
    --db-path ./literature_mining/knowledge_db.json

  # Search the local knowledge database
  python tools/literature_miner.py --query-db "diffusion" \\
    --db-path ./literature_mining/knowledge_db.json

  # Agent mode: ask questions against the knowledge database
  python tools/literature_miner.py --agent \\
    --question "What are key advances in protein diffusion?" \\
    --db-path ./literature_mining/knowledge_db.json

  # Agent mode: interactive Q&A session
  python tools/literature_miner.py --agent --agent-interactive \\
    --db-path ./literature_mining/knowledge_db.json

API Keys:
  - Tavily API:   set TAVILY_API_KEY env var or pass --tavily-api-key
  - DeepSeek API: set OPENAI_API_KEY env var or pass --deepseek-api-key
  - Or place a .env file in the repo root (auto-loaded):
      export TAVILY_API_KEY="your_key"
      export OPENAI_API_KEY="your_key"
      export OPENAI_BASE_URL="https://api.deepseek.com/v1"
"""

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import textwrap
import time
from datetime import datetime
from pathlib import Path
from urllib.parse import urlparse


# ---------------------------------------------------------------------------
# API defaults
# ---------------------------------------------------------------------------
TAVILY_API_URL = "https://api.tavily.com/search"
DEEPSEEK_BASE_URL = os.environ.get("OPENAI_BASE_URL", "https://api.deepseek.com/v1")
DEEPSEEK_MODEL = os.environ.get("OPENAI_MODEL", "deepseek-chat")

# ---------------------------------------------------------------------------
# Regex patterns for paper identification
# ---------------------------------------------------------------------------
DOI_PATTERN = re.compile(r"\b(10\.\d{4,9}/[-._;()/:A-Za-z0-9]+)\b")
ARXIV_ID_PATTERN = re.compile(r"arxiv\.org/abs/(\d+\.\d+(?:v\d+)?)")
ARXIV_RAW_PATTERN = re.compile(r"arxiv:\s*(\d+\.\d+(?:v\d+)?)", re.IGNORECASE)
PUBMED_PATTERN = re.compile(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)")
PMID_PATTERN = re.compile(r"PMID:\s*(\d+)", re.IGNORECASE)
PMCID_PATTERN = re.compile(r"(PMC\d+)", re.IGNORECASE)
PDF_URL_PATTERN = re.compile(r"\.pdf$", re.IGNORECASE)

# Keywords that suggest a URL is academic
ACADEMIC_DOMAINS = {
    "arxiv.org",
    "pubmed.ncbi.nlm.nih.gov",
    "ncbi.nlm.nih.gov",
    "pmc.ncbi.nlm.nih.gov",
    "nature.com",
    "science.org",
    "cell.com",
    "pnas.org",
    "academic.oup.com",
    "journals.plos.org",
    "pubs.acs.org",
    "onlinelibrary.wiley.com",
    "link.springer.com",
    "sciencedirect.com",
    "elifesciences.org",
    "biorxiv.org",
    "medrxiv.org",
    "chemrxiv.org",
    "researchgate.net",
    "semanticscholar.org",
    "dl.acm.org",
    "ieeexplore.ieee.org",
    ".edu",
    "mdpi.com",
    "frontiersin.org",
    "ncbi.nlm.nih.gov",
}


# ---------------------------------------------------------------------------
# .env file loading
# ---------------------------------------------------------------------------


def find_repo_root():
    """Find the repository root by looking for .git directory or .env file."""
    # Start from this script's directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Walk up until we find .git or .env
    current = script_dir
    for _ in range(5):
        if os.path.exists(os.path.join(current, ".git")) or os.path.exists(
            os.path.join(current, ".env")
        ):
            return current
        parent = os.path.dirname(current)
        if parent == current:
            break
        current = parent
    # Fallback to CWD
    return os.getcwd()


def load_env_file():
    """Load environment variables from a .env file (bash export format).

    Looks for .env in the repo root. Only sets variables not already
    present in the environment.
    """
    repo_root = find_repo_root()
    env_path = os.path.join(repo_root, ".env")
    if not os.path.exists(env_path):
        return

    env_re = re.compile(r'^(?:export\s+)?([A-Za-z_][A-Za-z0-9_]*)\s*=\s*"(.*)"\s*$')
    try:
        with open(env_path) as f:
            for line in f:
                m = env_re.match(line.strip())
                if m:
                    key, value = m.group(1), m.group(2)
                    if key not in os.environ:
                        os.environ[key] = value
    except Exception:
        pass  # Silently ignore .env parsing errors


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def get_api_key(arg_val, env_var):
    """Get API key from argument or environment variable."""
    if arg_val:
        return arg_val
    key = os.environ.get(env_var, "")
    if key:
        return key
    return None


def safe_import(module_name):
    """Lazy import with graceful fallback. Returns module or None."""
    try:
        return __import__(module_name)
    except ImportError:
        return None


def parse_args():
    parser = argparse.ArgumentParser(
        description="Literature Miner Tool - Search, download, and query academic papers",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples:
              # Search and download papers
              %(prog)s --query "protein diffusion models" --max-downloads 5

              # Search the local knowledge database
              %(prog)s --query-db "diffusion" --db-path ./literature_mining/knowledge_db.json

              # Agent mode: ask questions against the knowledge database
              %(prog)s --agent --question "What are key advances in protein diffusion?"

              # Interactive agent session
              %(prog)s --agent --agent-interactive --db-path ./literature_mining/knowledge_db.json
        """),
    )

    # --- Mode selection (mutually exclusive) ---
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument(
        "--query", "-q", help="Topic or sentence to search the web for"
    )
    mode_group.add_argument(
        "--query-db", help="Search the local knowledge database for a keyword"
    )
    mode_group.add_argument(
        "--agent",
        action="store_true",
        help="Enable agent/RAG mode for Q&A over the knowledge database",
    )
    mode_group.add_argument(
        "--pdf-dir",
        help="Path to a directory of PDF files to ingest into the knowledge database",
    )

    # --- Web search options ---
    parser.add_argument(
        "--output-dir",
        "-o",
        default="./literature_mining",
        help="Output directory (default: ./literature_mining)",
    )
    parser.add_argument(
        "--tavily-api-key", help="Tavily API key (or set TAVILY_API_KEY env var)"
    )
    parser.add_argument(
        "--deepseek-api-key", help="DeepSeek API key (or set OPENAI_API_KEY env var)"
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=20,
        help="Max search results to process (default: 20)",
    )
    parser.add_argument(
        "--max-downloads",
        type=int,
        default=10,
        help="Max papers to download (default: 10)",
    )
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="Skip PDF download, only collect metadata",
    )
    parser.add_argument(
        "--no-summarize", action="store_true", help="Skip LLM summarization of papers"
    )
    parser.add_argument(
        "--scihub-mirror",
        default="https://sci-hub.se",
        help="Sci-Hub mirror URL (default: https://sci-hub.se)",
    )

    # --- Agent mode options ---
    parser.add_argument(
        "--question", help="Question to ask in agent mode (used with --agent)"
    )
    parser.add_argument(
        "--agent-top-k",
        type=int,
        default=5,
        help="Top papers for context in agent mode (default: 5)",
    )
    parser.add_argument(
        "--agent-interactive",
        action="store_true",
        help="Interactive Q&A session (used with --agent)",
    )

    # --- Database options ---
    parser.add_argument(
        "--db-path", default=None, help="Path to knowledge database JSON file"
    )

    # --- PDF directory options ---
    parser.add_argument(
        "--pdf-max-chars",
        type=int,
        default=8000,
        help="Max characters to extract from each PDF (default: 8000)",
    )
    parser.add_argument(
        "--no-pdf-skip-existing",
        action="store_false",
        dest="pdf_skip_existing",
        help="Re-process PDFs even if already in the knowledge database",
    )

    # --- Output options ---
    parser.add_argument(
        "--output-format",
        default="table",
        choices=["table", "json", "csv"],
        help="Output format (default: table)",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Web Search (Tavily API)
# ---------------------------------------------------------------------------


def search_tavily(query, api_key, max_results=20):
    """Search Tavily API and return a list of result dicts."""
    try:
        import requests
    except ImportError as e:
        print(
            f"Error: 'requests' package is required ({e}). Install with: pip install requests",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Searching Tavily for: '{query}'...")
    try:
        resp = requests.post(
            TAVILY_API_URL,
            json={
                "query": query,
                "max_results": max_results,
                "search_depth": "advanced",
                "include_answer": False,
                "include_raw_content": False,
            },
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {api_key}",
            },
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json()
        results = data.get("results", [])
        print(f"  Found {len(results)} results from Tavily.")
        return results
    except Exception as e:
        print(f"Error searching Tavily: {e}", file=sys.stderr)
        return []


# ---------------------------------------------------------------------------
# Paper identification
# ---------------------------------------------------------------------------


def is_academic_url(url):
    """Check if a URL is likely to be an academic source."""
    try:
        parsed = urlparse(url)
        domain = parsed.netloc.lower()
        for ad in ACADEMIC_DOMAINS:
            if ad in domain:
                return True
        # Also check if the URL contains academic path patterns
        path = parsed.path.lower()
        if any(
            kw in path
            for kw in [
                "/abs/",
                "/doi/",
                "/article/",
                "/paper/",
                "/publication/",
                "/pmc/",
                "/pubmed/",
            ]
        ):
            return True
        return False
    except Exception:
        return False


def extract_identifiers(url, content_snippet=""):
    """Extract paper identifiers (DOI, arXiv ID, PubMed ID, PMCID) from URL and snippet.

    Returns a dict with keys: doi, arxiv_id, pubmed_id, pmcid.
    """
    text = url + " " + (content_snippet or "")
    ids = {"doi": None, "arxiv_id": None, "pubmed_id": None, "pmcid": None}

    m = DOI_PATTERN.search(text)
    if m:
        ids["doi"] = m.group(1)

    m = ARXIV_ID_PATTERN.search(text)
    if m:
        ids["arxiv_id"] = m.group(1)
    else:
        m = ARXIV_RAW_PATTERN.search(text)
        if m:
            ids["arxiv_id"] = m.group(1)

    m = PUBMED_PATTERN.search(text)
    if m:
        ids["pubmed_id"] = m.group(1)
    else:
        m = PMID_PATTERN.search(text)
        if m:
            ids["pubmed_id"] = m.group(1)

    m = PMCID_PATTERN.search(text)
    if m:
        ids["pmcid"] = m.group(1)

    return ids


# ---------------------------------------------------------------------------
# PDF Download
# ---------------------------------------------------------------------------


def find_pdf_link_from_page(url):
    """Fetch an article page and look for the PDF link in meta tags or common patterns.

    Returns the PDF URL if found, None otherwise.
    """
    requests = safe_import("requests")
    if not requests:
        return None
    try:
        resp = requests.get(
            url,
            timeout=30,
            allow_redirects=True,
            headers={"User-Agent": "Mozilla/5.0 (compatible; LiteratureMiner/1.0)"},
        )
        html = resp.text
        # Look for citation_pdf_url meta tag
        m = re.search(
            r'<meta\s+name="citation_pdf_url"\s+content="([^"]+)"', html, re.IGNORECASE
        )
        if m:
            return m.group(1)
        # Look for PDF links in page
        m = re.search(r'href="([^"]+\.pdf)"', html, re.IGNORECASE)
        if m:
            pdf_url = m.group(1)
            if pdf_url.startswith("/"):
                base = "{0.scheme}://{0.netloc}".format(urlparse(resp.url))
                pdf_url = base + pdf_url
            return pdf_url
        return None
    except Exception:
        return None


def build_download_urls(identifiers, source_url, scihub_mirror="https://sci-hub.se"):
    """Build a prioritized list of download URLs for a paper.

    Returns list of (source_label, url) tuples.
    """
    urls = []

    if identifiers.get("arxiv_id"):
        arxiv_id = identifiers["arxiv_id"]
        urls.append(("arXiv", f"https://arxiv.org/pdf/{arxiv_id}.pdf"))

    if identifiers.get("pmcid"):
        pmcid = identifiers["pmcid"]
        # Europe PMC first (most reliable for open-access PDFs)
        urls.append(
            (
                "Europe PMC",
                f"https://www.ebi.ac.uk/europepmc/webservices/rest/pmc/{pmcid}/fullTextPDF",
            )
        )
        urls.append(("Europe PMC (alt)", f"https://europepmc.org/articles/{pmcid}/pdf"))
        # NCBI PMC as fallback
        urls.append(
            (
                "PubMed Central",
                f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/main.pdf",
            )
        )

    if identifiers.get("pubmed_id"):
        pubmed_id = identifiers["pubmed_id"]
        # Europe PMC supports PubMed IDs too
        urls.append(
            (
                "Europe PMC (PMID)",
                f"https://www.ebi.ac.uk/europepmc/webservices/rest/MED/{pubmed_id}/fullTextPDF",
            )
        )
        urls.append(
            (
                "Europe PMC (PMID alt)",
                f"https://europepmc.org/article/MED/{pubmed_id}?pdf=render",
            )
        )

    if identifiers.get("doi"):
        doi = identifiers["doi"]
        scihub = scihub_mirror.rstrip("/")
        # Try DOI resolution via doi.org (some publishers serve PDFs directly)
        urls.append(("DOI (doi.org)", f"https://doi.org/{doi}"))
        # Multiple Sci-Hub mirrors
        urls.append(("Sci-Hub", f"{scihub}/{doi}"))
        urls.append(("Sci-Hub (ru)", f"https://sci-hub.ru/{doi}"))
        urls.append(("Sci-Hub (st)", f"https://sci-hub.st/{doi}"))
        # For biorxiv DOIs, try direct PDF
        if "biorxiv" in source_url.lower() or "/10.1101/" in doi:
            urls.append(("bioRxiv", f"https://www.biorxiv.org/content/{doi}.full.pdf"))

    # Direct URL if it points to a PDF
    if PDF_URL_PATTERN.search(source_url):
        urls.append(("Direct PDF", source_url))

    # For academic domains without explicit IDs, try source URL anyway
    if (
        not identifiers.get("doi")
        and not identifiers.get("arxiv_id")
        and not identifiers.get("pmcid")
    ):
        if is_academic_url(source_url) and not any(u for _, u in urls):
            urls.append(("Source page", source_url))

    return urls


def download_pdf(url, dest_path, label):
    """Download a PDF from a URL. Returns True on success."""
    requests = safe_import("requests")
    if not requests:
        return False

    try:
        print(f"    Trying {label}: {url[:100]}...")
        resp = requests.get(
            url,
            timeout=60,
            allow_redirects=True,
            headers={"User-Agent": "Mozilla/5.0 (compatible; LiteratureMiner/1.0)"},
        )
        if resp.status_code == 200:
            content_type = resp.headers.get("Content-Type", "")
            content = resp.content
            # Check magic bytes first (most reliable), then Content-Type, then URL extension
            is_pdf_magic = content[:5] == b"%PDF-"
            is_pdf_mime = "application/pdf" in content_type
            is_pdf_url = url.endswith(".pdf")
            if is_pdf_magic:
                if len(content) < 1000:
                    print(f"    PDF too small ({len(content)} bytes), skipping")
                    return False
                with open(dest_path, "wb") as f:
                    f.write(content)
                print(f"    Downloaded: {dest_path} ({len(content)} bytes)")
                return True
            elif is_pdf_mime or is_pdf_url:
                # Content-Type says PDF or URL ends in .pdf, but magic bytes mismatch
                print(
                    f"    Server claimed PDF but content is not (first bytes: {content[:20]!r})"
                )
                return False
            else:
                print(f"    Not a PDF (Content-Type: {content_type})")
                return False
        else:
            print(f"    HTTP {resp.status_code}")
            return False
    except Exception as e:
        print(f"    Error: {e}")
        return False


def download_paper(identifiers, source_url, paper_dir, scihub_mirror):
    """Download a single paper's PDF.

    Returns the local PDF path on success, None on failure.
    """
    urls = build_download_urls(identifiers, source_url, scihub_mirror)
    if not urls:
        print("    No download URLs identified")
        return None

    pdf_name = None
    if identifiers.get("doi"):
        pdf_name = identifiers["doi"].replace("/", "_") + ".pdf"
    elif identifiers.get("arxiv_id"):
        pdf_name = f"arxiv_{identifiers['arxiv_id']}.pdf"
    elif identifiers.get("pubmed_id"):
        pdf_name = f"pubmed_{identifiers['pubmed_id']}.pdf"
    else:
        pdf_name = hashlib.md5(source_url.encode()).hexdigest()[:12] + ".pdf"

    dest_path = os.path.join(paper_dir, pdf_name)
    if os.path.exists(dest_path):
        print(f"    PDF already exists: {dest_path}")
        return dest_path

    for label, url in urls:
        if download_pdf(url, dest_path, label):
            return dest_path

    # Fallback: scrape the article page for PDF link
    if source_url:
        print(f"    Scraping article page for PDF link: {source_url[:80]}...")
        pdf_url = find_pdf_link_from_page(source_url)
        if pdf_url:
            print(f"    Found PDF link: {pdf_url[:100]}...")
            if download_pdf(pdf_url, dest_path, "Page-scraped PDF"):
                return dest_path

    return None


# ---------------------------------------------------------------------------
# PDF text extraction
# ---------------------------------------------------------------------------


def is_valid_pdf(path):
    """Quick check if a file starts with the PDF magic bytes."""
    try:
        with open(path, "rb") as f:
            return f.read(5) == b"%PDF-"
    except Exception:
        return False


def extract_pdf_text(pdf_path, max_chars=8000):
    """Extract text from a PDF file. Returns first max_chars characters.

    Returns None if the file is not a valid PDF or has no extractable text.
    """
    if not is_valid_pdf(pdf_path):
        return None

    import contextlib

    # Suppress pypdf's stderr noise
    stderr_null = open(os.devnull, "w")

    def _try_extract(reader):
        text_parts = []
        total = 0
        for page in reader.pages:
            page_text = page.extract_text()
            if page_text:
                text_parts.append(page_text)
                total += len(page_text)
                if total >= max_chars:
                    break
        return "\n".join(text_parts)[:max_chars] if text_parts else None

    # Try PyPDF2 first
    PyPDF2 = safe_import("PyPDF2")
    if PyPDF2:
        try:
            with contextlib.redirect_stderr(stderr_null):
                with contextlib.redirect_stdout(stderr_null):
                    reader = PyPDF2.PdfReader(pdf_path)
                    result = _try_extract(reader)
            if result:
                return result
        except Exception:
            pass

    # Try pdfplumber as fallback
    pdfplumber = safe_import("pdfplumber")
    if pdfplumber:
        try:
            with contextlib.redirect_stderr(stderr_null):
                with contextlib.redirect_stdout(stderr_null):
                    with pdfplumber.open(pdf_path) as pdf:
                        text_parts = []
                        total = 0
                        for page in pdf.pages:
                            page_text = page.extract_text()
                            if page_text:
                                text_parts.append(page_text)
                                total += len(page_text)
                                if total >= max_chars:
                                    break
                        result = (
                            "\n".join(text_parts)[:max_chars] if text_parts else None
                        )
            if result:
                return result
        except Exception:
            pass

    # Try pypdf (modern PyPDF2 fork) as last fallback
    pypdf = safe_import("pypdf")
    if pypdf:
        try:
            with contextlib.redirect_stderr(stderr_null):
                with contextlib.redirect_stdout(stderr_null):
                    reader = pypdf.PdfReader(pdf_path)
                    result = _try_extract(reader)
            if result:
                return result
        except Exception:
            pass

    # Cleanup
    try:
        stderr_null.close()
    except Exception:
        pass

    return None


# ---------------------------------------------------------------------------
# LLM (DeepSeek) helpers
# ---------------------------------------------------------------------------


def call_deepseek(prompt, api_key, system_prompt=None, temperature=0.3):
    """Call the DeepSeek API (OpenAI-compatible) and return the response text."""
    try:
        import requests
    except ImportError as e:
        print(f"Error: 'requests' package is required ({e}).", file=sys.stderr)
        return None

    messages = []
    if system_prompt:
        messages.append({"role": "system", "content": system_prompt})
    messages.append({"role": "user", "content": prompt})

    url = f"{DEEPSEEK_BASE_URL}/chat/completions"
    try:
        resp = requests.post(
            url,
            json={
                "model": DEEPSEEK_MODEL,
                "messages": messages,
                "temperature": temperature,
                "max_tokens": 2000,
            },
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {api_key}",
            },
            timeout=120,
        )
        resp.raise_for_status()
        data = resp.json()
        return data["choices"][0]["message"]["content"].strip()
    except Exception as e:
        print(f"  DeepSeek API error: {e}", file=sys.stderr)
        return None


def summarize_paper(title, source_url, abstract, pdf_text, api_key):
    """Use DeepSeek to summarize a paper and extract key findings."""
    text_content = ""
    if abstract:
        text_content += f"Abstract: {abstract}\n\n"
    if pdf_text:
        text_content += f"Full text excerpt: {pdf_text[:5000]}"

    if not text_content.strip():
        return None, None, []

    system = (
        "You are a scientific literature analyst. "
        "Provide concise, accurate summaries in JSON format."
    )

    prompt = f"""Analyze this scientific paper and provide:
1. A 2-3 sentence summary of the main contribution
2. Up to 5 key findings as a list
3. Up to 5 relevant keywords

Title: {title}
Source: {source_url}

{text_content}

Respond in this exact JSON format:
{{
  "summary": "2-3 sentence summary here",
  "key_findings": ["finding 1", "finding 2", ...],
  "keywords": ["keyword1", "keyword2", ...]
}}"""

    response = call_deepseek(prompt, api_key, system_prompt=system)
    if not response:
        return None, None, []

    # Try to parse JSON from response
    try:
        # Find JSON block in response
        json_match = re.search(r"\{[\s\S]*\}", response)
        if json_match:
            data = json.loads(json_match.group())
            return (
                data.get("summary", ""),
                data.get("key_findings", []),
                data.get("keywords", []),
            )
    except (json.JSONDecodeError, KeyError):
        pass

    # Fallback: use the raw response as summary
    return response, [], []


# ---------------------------------------------------------------------------
# Knowledge Database
# ---------------------------------------------------------------------------


def load_knowledge_db(db_path):
    """Load the knowledge database from a JSON file. Returns dict or None."""
    if not os.path.exists(db_path):
        return None
    try:
        with open(db_path) as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError) as e:
        print(f"Error loading knowledge database: {e}", file=sys.stderr)
        return None


def save_knowledge_db(db, db_path):
    """Save the knowledge database to a JSON file."""
    os.makedirs(os.path.dirname(db_path) or ".", exist_ok=True)
    with open(db_path, "w") as f:
        json.dump(db, f, indent=2, ensure_ascii=False)


def init_knowledge_db():
    """Create a fresh knowledge database structure."""
    return {
        "papers": [],
        "search_history": [],
    }


def add_paper_to_db(db, paper_entry):
    """Add a paper entry to the database. Avoid duplicates by DOI/arXiv ID."""
    existing_ids = set()
    for p in db["papers"]:
        if p.get("doi"):
            existing_ids.add(("doi", p["doi"]))
        if p.get("arxiv_id"):
            existing_ids.add(("arxiv", p["arxiv_id"]))
        if p.get("pubmed_id"):
            existing_ids.add(("pubmed", p["pubmed_id"]))

    # Check duplicates
    if paper_entry.get("doi") and ("doi", paper_entry["doi"]) in existing_ids:
        return False
    if (
        paper_entry.get("arxiv_id")
        and ("arxiv", paper_entry["arxiv_id"]) in existing_ids
    ):
        return False
    if (
        paper_entry.get("pubmed_id")
        and ("pubmed", paper_entry["pubmed_id"]) in existing_ids
    ):
        return False

    db["papers"].append(paper_entry)
    return True


def search_db(db, query, top_k=None):
    """Search the knowledge database for papers matching a query.

    Searches across: title, abstract, summary, key_findings, keywords, authors.
    Returns a list of (paper, score) tuples sorted by relevance.
    """
    if not db or not db.get("papers"):
        return []

    query_lower = query.lower()
    query_terms = query_lower.split()

    scored = []
    for paper in db["papers"]:
        score = 0

        # Build searchable text with field weights
        search_fields = [
            (paper.get("title") or "", 10),  # title: weight 10
            (paper.get("keywords") or [], 8),  # keywords: weight 8
            (paper.get("summary") or "", 6),  # summary: weight 6
            (paper.get("key_findings") or [], 5),  # findings: weight 5
            (paper.get("abstract") or "", 3),  # abstract: weight 3
            (paper.get("search_query") or "", 2),  # original query: weight 2
        ]

        for field, weight in search_fields:
            if isinstance(field, list):
                field_text = " ".join(field).lower()
            else:
                field_text = field.lower()

            # Exact phrase match (higher boost)
            if query_lower in field_text:
                score += weight * 5
            # Individual term matches
            for term in query_terms:
                if term in field_text:
                    score += weight

        if score > 0:
            scored.append((paper, score))

    scored.sort(key=lambda x: x[1], reverse=True)
    if top_k:
        scored = scored[:top_k]
    return scored


# ---------------------------------------------------------------------------
# Agent Mode (RAG Q&A)
# ---------------------------------------------------------------------------


def build_agent_context(papers_with_scores):
    """Build a context prompt from retrieved papers."""
    context_parts = []
    for i, (paper, score) in enumerate(papers_with_scores, 1):
        title = paper.get("title", "Unknown Title")
        authors = ", ".join(paper.get("authors", [])) or "Unknown"
        summary = paper.get("summary", "No summary available")
        findings = paper.get("key_findings", [])
        doi = paper.get("doi", "")
        arxiv = paper.get("arxiv_id", "")
        source = paper.get("source_url", "")

        part = f"[Paper {i}]\n"
        part += f"Title: {title}\n"
        part += f"Authors: {authors}\n"
        if doi:
            part += f"DOI: https://doi.org/{doi}\n"
        if arxiv:
            part += f"arXiv: https://arxiv.org/abs/{arxiv}\n"
        if source:
            part += f"Source: {source}\n"
        part += f"Summary: {summary}\n"
        if findings:
            part += "Key Findings:\n"
            for f in findings:
                part += f"  - {f}\n"

        context_parts.append(part)

    return "\n".join(context_parts)


def agent_query(db, question, api_key, top_k=5):
    """Run a RAG query: retrieve relevant papers and generate an answer.

    Returns (answer, retrieved_papers) tuple.
    """
    if not db or not db.get("papers"):
        return (
            "No papers in the knowledge database. Run a search first with --query.",
            [],
        )

    # Retrieve top-K relevant papers
    results = search_db(db, question, top_k=top_k)
    if not results:
        return ("No relevant papers found in the database for this question.", [])

    print(f"  Retrieved {len(results)} relevant papers for context")

    # Build context
    context = build_agent_context(results)

    # Build the agent prompt
    system = (
        "You are a research assistant with access to a literature database. "
        "Answer questions based on the provided paper summaries. "
        "Be concise and cite specific papers using [Paper N] notation. "
        "If the papers do not contain enough information to answer, say so."
    )

    prompt = f"""Based on the following paper summaries, answer the question.

{context}

QUESTION: {question}

Answer concisely, citing specific papers with [Paper N]. Include key details from the papers."""

    answer = call_deepseek(prompt, api_key, system_prompt=system, temperature=0.3)
    if not answer:
        return (
            "Failed to generate answer from DeepSeek API.",
            [(p, s) for p, s in results],
        )

    return answer, [(p, s) for p, s in results]


def interactive_agent(db, api_key, top_k):
    """Run an interactive agent Q&A session."""
    if not db or not db.get("papers"):
        print("No papers in the knowledge database. Run a search first with --query.")
        return

    n_papers = len(db["papers"])
    print(f"\n{'=' * 60}")
    print(f"INTERACTIVE AGENT MODE")
    print(f"{'=' * 60}")
    print(f"Knowledge database: {n_papers} papers loaded")
    print(f"Context window: top {top_k} papers per question")
    print(f"Type 'exit' or 'quit' to end the session, 'help' for commands.")
    print()

    while True:
        try:
            question = input("Ask> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nExiting...")
            break

        if not question:
            continue
        if question.lower() in ("exit", "quit"):
            print("Goodbye.")
            break
        if question.lower() == "help":
            print("  Commands:")
            print("    exit/quit  - End the session")
            print("    stats      - Show database statistics")
            print("    papers     - List paper titles in the database")
            continue
        if question.lower() == "stats":
            print(f"  Total papers: {n_papers}")
            keywords = set()
            for p in db["papers"]:
                for kw in p.get("keywords", []):
                    keywords.add(kw)
            print(f"  Unique keywords: {len(keywords)}")
            continue
        if question.lower() == "papers":
            for i, p in enumerate(db["papers"], 1):
                print(f"  [{i}] {p.get('title', 'Untitled')[:100]}")
            continue

        print(f"  Searching database...")
        answer, retrieved = agent_query(db, question, api_key, top_k=top_k)
        print(f"\n{answer}\n")
        if retrieved:
            print("Sources:")
            for i, (paper, score) in enumerate(retrieved, 1):
                title = paper.get("title", "Unknown")[:80]
                doi = paper.get("doi", "")
                src = f"https://doi.org/{doi}" if doi else paper.get("source_url", "")
                print(f"  [Paper {i}] {title}")
                if src:
                    print(f"           {src}")
            print()


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------


def print_results_table(papers):
    """Print papers as a formatted table."""
    if not papers:
        print("No results to display.")
        return

    print(f"\n{'=' * 100}")
    print(f"RESULTS: {len(papers)} paper(s)")
    print(f"{'=' * 100}")

    for i, (paper, score) in enumerate(papers, 1):
        title = paper.get("title", "Untitled")[:90]
        doi = paper.get("doi", "")
        arxiv = paper.get("arxiv_id", "")
        source = paper.get("source_url", "")
        pdf = paper.get("pdf_path", "")

        print(f"\n--- Paper {i} (score: {score}) ---")
        print(f"  Title:    {title}")
        if doi:
            print(f"  DOI:      https://doi.org/{doi}")
        if arxiv:
            print(f"  arXiv:    https://arxiv.org/abs/{arxiv}")
        if source:
            print(f"  Source:   {source}")
        if pdf:
            print(f"  PDF:      {pdf}")
        abstract = paper.get("abstract", "")
        if abstract:
            print(f"  Abstract: {abstract[:200]}{'...' if len(abstract) > 200 else ''}")
        summary = paper.get("summary", "")
        if summary:
            print(f"  Summary:  {summary[:200]}{'...' if len(summary) > 200 else ''}")
        findings = paper.get("key_findings", [])
        if findings:
            print(f"  Findings: ")
            for f in findings[:3]:
                print(f"    - {f}")


def print_papers_json(papers):
    """Print papers as JSON."""
    output = [p for p, _ in papers]
    print(json.dumps(output, indent=2, ensure_ascii=False))


def print_papers_csv(papers):
    """Print papers as CSV."""
    if not papers:
        return
    fieldnames = [
        "title",
        "doi",
        "arxiv_id",
        "pubmed_id",
        "source_url",
        "pdf_path",
        "abstract",
        "summary",
        "keywords",
    ]
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    for paper, _ in papers:
        row = {k: paper.get(k, "") for k in fieldnames}
        if isinstance(row.get("keywords"), list):
            row["keywords"] = "; ".join(row["keywords"])
        writer.writerow(row)


# ---------------------------------------------------------------------------
# Mode A: Web Search & Download
# ---------------------------------------------------------------------------


def mode_search_and_download(args):
    """Run the web search, download, and summarize pipeline."""
    # --- API key checks ---
    tavily_key = get_api_key(args.tavily_api_key, "TAVILY_API_KEY")
    if not tavily_key:
        print(
            "Error: Tavily API key required. Set TAVILY_API_KEY env var or use --tavily-api-key.",
            file=sys.stderr,
        )
        sys.exit(1)

    deepseek_key = None
    if not args.no_summarize:
        deepseek_key = get_api_key(args.deepseek_api_key, "OPENAI_API_KEY")
        if not deepseek_key:
            print(
                "Warning: No DeepSeek API key found. Summarization disabled. "
                "Set OPENAI_API_KEY env var or use --deepseek-api-key."
            )
            args.no_summarize = True

    # --- Setup ---
    output_dir = os.path.abspath(args.output_dir)
    paper_dir = os.path.join(output_dir, "papers")
    db_path = args.db_path or os.path.join(output_dir, "knowledge_db.json")

    os.makedirs(output_dir, exist_ok=True)
    if not args.no_download:
        os.makedirs(paper_dir, exist_ok=True)

    print(f"{'=' * 60}")
    print(f"LITERATURE MINER - Web Search & Download")
    print(f"{'=' * 60}")
    print(f"Query:        {args.query}")
    print(f"Output dir:   {output_dir}")
    print(f"Max results:  {args.max_results}")
    print(f"Max download: {args.max_downloads}")
    print(f"Download:     {not args.no_download}")
    print(f"Summarize:    {not args.no_summarize}")

    # --- Phase 1: Search (multiple queries targeting OA sources) ---
    print(f"\n{'=' * 60}")
    print("PHASE 1: Web Search")
    print(f"{'=' * 60}")

    # Build queries: bare query + site-specific OA queries
    queries = [args.query]
    oa_sites = ["arxiv.org", "biorxiv.org", "frontiersin.org", "mdpi.com", "plos.org"]
    for site in oa_sites:
        queries.append(f"{args.query} site:{site}")

    all_results = []
    seen_urls = set()
    for q in queries:
        results = search_tavily(q, tavily_key, max(args.max_results // 3, 10))
        for r in results:
            url = r.get("url", "")
            if url and url not in seen_urls:
                seen_urls.add(url)
                all_results.append(r)
        print(f"  Query '{q[:60]}...' → {len(results)} results")

    results = all_results[: args.max_results]
    print(f"  Total unique results: {len(results)}")

    if not results:
        print("No search results found.")
        return

    # --- Phase 2: Identify papers ---
    print(f"\n{'=' * 60}")
    print("PHASE 2: Paper Identification")
    print(f"{'=' * 60}")

    papers_meta = []
    for r in results:
        url = r.get("url", "")
        title = r.get("title", "Untitled")
        content = r.get("content", "")
        score = r.get("score", 0)

        if not url:
            continue

        identifiers = extract_identifiers(url, content)
        has_identifier = any(identifiers.values())
        is_academic = is_academic_url(url)

        if has_identifier or is_academic:
            paper = {
                "id": hashlib.md5(url.encode()).hexdigest()[:12],
                "title": title,
                "source_url": url,
                "abstract": content[:1000] if content else "",
                "score": score,
                "doi": identifiers["doi"],
                "arxiv_id": identifiers["arxiv_id"],
                "pubmed_id": identifiers["pubmed_id"],
                "pmcid": identifiers["pmcid"],
            }
            papers_meta.append(paper)
            id_str = (
                identifiers["doi"]
                or identifiers["arxiv_id"]
                or identifiers["pubmed_id"]
                or "no-id"
            )
            print(f"  [{len(papers_meta)}] {title[:80]}...")
            print(f"      ID: {id_str} | URL: {url[:80]}...")

    print(f"\nIdentified {len(papers_meta)} academic papers.")

    if not papers_meta:
        return

    # --- Phase 3: Download ---
    pdf_map = {}  # paper_id -> pdf_path
    if not args.no_download:
        print(f"\n{'=' * 60}")
        print(f"PHASE 3: PDF Download (max: {args.max_downloads})")
        print(f"{'=' * 60}")

        downloads = 0
        for paper in papers_meta:
            if downloads >= args.max_downloads:
                print(f"  Reached download limit ({args.max_downloads}).")
                break

            identifiers = {
                "doi": paper["doi"],
                "arxiv_id": paper["arxiv_id"],
                "pubmed_id": paper["pubmed_id"],
                "pmcid": paper["pmcid"],
            }

            print(f"  [{downloads + 1}] {paper['title'][:80]}...")
            pdf_path = download_paper(
                identifiers, paper["source_url"], paper_dir, args.scihub_mirror
            )
            if pdf_path:
                pdf_map[paper["id"]] = pdf_path
                downloads += 1

        print(f"\nDownloaded {downloads} PDF(s) to {paper_dir}/")

    # --- Phase 4: Summarize ---
    if not args.no_summarize and deepseek_key:
        print(f"\n{'=' * 60}")
        print("PHASE 4: LLM Summarization")
        print(f"{'=' * 60}")

        for paper in papers_meta:
            print(f"  Summarizing: {paper['title'][:80]}...")
            pdf_path = pdf_map.get(paper["id"])
            pdf_text = None
            if pdf_path:
                print(f"    Extracting PDF text...")
                pdf_text = extract_pdf_text(pdf_path)

            summary, findings, keywords = summarize_paper(
                paper["title"],
                paper["source_url"],
                paper.get("abstract", ""),
                pdf_text,
                deepseek_key,
            )
            paper["summary"] = summary or ""
            paper["key_findings"] = findings
            paper["keywords"] = keywords

            if summary:
                print(f"    Summary: {summary[:120]}...")
            if findings:
                print(f"    Findings: {len(findings)}")
            if keywords:
                print(f"    Keywords: {', '.join(keywords[:5])}")

    # --- Phase 5: Store ---
    print(f"\n{'=' * 60}")
    print("PHASE 5: Store to Knowledge Database")
    print(f"{'=' * 60}")

    db = load_knowledge_db(db_path)
    if db is None:
        db = init_knowledge_db()

    new_count = 0
    for paper in papers_meta:
        pdf_path = pdf_map.get(paper["id"])
        entry = {
            "id": paper["id"],
            "title": paper["title"],
            "authors": paper.get("authors", []),
            "abstract": paper.get("abstract", ""),
            "source_url": paper["source_url"],
            "doi": paper["doi"],
            "arxiv_id": paper["arxiv_id"],
            "pubmed_id": paper["pubmed_id"],
            "pdf_path": pdf_path,
            "keywords": paper.get("keywords", []),
            "summary": paper.get("summary", ""),
            "key_findings": paper.get("key_findings", []),
            "search_query": args.query,
            "date_added": datetime.now().strftime("%Y-%m-%d"),
        }
        if add_paper_to_db(db, entry):
            new_count += 1

    # Record search history
    db["search_history"].append(
        {
            "query": args.query,
            "date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "results_count": len(papers_meta),
            "new_papers": new_count,
        }
    )

    save_knowledge_db(db, db_path)
    print(f"  Added {new_count} new papers to database")
    print(f"  Database: {db_path} ({len(db['papers'])} total papers)")

    # --- Phase 6: Report ---
    print(f"\n{'=' * 60}")
    print("PHASE 6: Report")
    print(f"{'=' * 60}")

    if args.output_format == "json":
        print_papers_json([(p, p.get("score", 0)) for p in papers_meta])
    elif args.output_format == "csv":
        print_papers_csv([(p, p.get("score", 0)) for p in papers_meta])
    else:
        print_results_table([(p, p.get("score", 0)) for p in papers_meta])

    print(f"\nLiterature mining complete.")
    print(f"  Papers identified: {len(papers_meta)}")
    print(f"  PDFs downloaded:   {len(pdf_map)}")
    print(f"  Database:          {db_path}")
    print(f"\nTo query the database:")
    print(
        f'  python tools/literature_miner.py --query-db "<keyword>" --db-path {db_path}'
    )
    print(f"To ask questions in agent mode:")
    print(
        f'  python tools/literature_miner.py --agent --question "<question>" --db-path {db_path}'
    )


# ---------------------------------------------------------------------------
# Mode B: Database Search
# ---------------------------------------------------------------------------


def mode_query_db(args):
    """Search the local knowledge database."""
    db_path = args.db_path
    if not db_path:
        print("Error: --db-path is required for database search.", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(db_path):
        print(f"Error: Database not found at {db_path}", file=sys.stderr)
        sys.exit(1)

    db = load_knowledge_db(db_path)
    if not db or not db.get("papers"):
        print("Database is empty. Run a search first with --query.")
        return

    query = args.query_db
    print(f"\nSearching database for: '{query}'")
    print(f"Database: {db_path} ({len(db['papers'])} papers)")

    results = search_db(db, query)
    if not results:
        print(f"\nNo papers found matching '{query}'.")
        return

    if args.output_format == "json":
        print_papers_json(results)
    elif args.output_format == "csv":
        print_papers_csv(results)
    else:
        print_results_table(results)


# ---------------------------------------------------------------------------
# Mode C: Agent RAG Q&A
# ---------------------------------------------------------------------------


def mode_agent(args):
    """Run agent/RAG mode for Q&A over the knowledge database."""
    db_path = args.db_path
    if not db_path:
        print("Error: --db-path is required for agent mode.", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(db_path):
        print(f"Error: Database not found at {db_path}", file=sys.stderr)
        print("Run a search first with --query to build the knowledge database.")
        sys.exit(1)

    db = load_knowledge_db(db_path)
    if not db or not db.get("papers"):
        print("Database is empty. Run a search first with --query.")
        return

    deepseek_key = get_api_key(args.deepseek_api_key, "OPENAI_API_KEY")
    if not deepseek_key:
        print(
            "Error: DeepSeek API key required for agent mode. "
            "Set OPENAI_API_KEY env var or use --deepseek-api-key.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.agent_interactive:
        interactive_agent(db, deepseek_key, args.agent_top_k)
        return

    if not args.question:
        print(
            "Error: --question is required for agent mode (unless using --agent-interactive).",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"\n{'=' * 60}")
    print(f"AGENT MODE - RAG Q&A")
    print(f"{'=' * 60}")
    print(f"Database: {db_path} ({len(db['papers'])} papers)")
    print(f"Question: {args.question}")
    print(f"Context:  top {args.agent_top_k} papers")

    answer, retrieved = agent_query(db, args.question, deepseek_key, args.agent_top_k)

    print(f"\n{'=' * 60}")
    print("ANSWER")
    print(f"{'=' * 60}")
    print(f"\n{answer}\n")

    if retrieved:
        print(f"{'=' * 60}")
        print("SOURCES")
        print(f"{'=' * 60}")
        for i, (paper, score) in enumerate(retrieved, 1):
            title = paper.get("title", "Unknown")
            doi = paper.get("doi", "")
            arxiv = paper.get("arxiv_id", "")
            src = paper.get("source_url", "")

            print(f"\n[Paper {i}] (relevance: {score})")
            print(f"  Title: {title}")
            if doi:
                print(f"  DOI:   https://doi.org/{doi}")
            if arxiv:
                print(f"  arXiv: https://arxiv.org/abs/{arxiv}")
            if src:
                print(f"  URL:   {src}")


# ---------------------------------------------------------------------------
# Mode D: PDF Directory Ingestion
# ---------------------------------------------------------------------------


def hash_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(65536)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def mode_pdf_dir(args):
    pdf_dir = os.path.abspath(args.pdf_dir)
    if not os.path.isdir(pdf_dir):
        print(f"Error: PDF directory not found: {pdf_dir}", file=sys.stderr)
        sys.exit(1)

    deepseek_key = get_api_key(args.deepseek_api_key, "OPENAI_API_KEY")
    if not deepseek_key:
        print(
            "Error: DeepSeek API key required for PDF ingestion. "
            "Set OPENAI_API_KEY env var or use --deepseek-api-key.",
            file=sys.stderr,
        )
        sys.exit(1)

    db_path = args.db_path or os.path.join(
        os.path.dirname(pdf_dir), "knowledge_db.json"
    )

    db = load_knowledge_db(db_path)
    if db is None:
        db = init_knowledge_db()

    existing_hashes = set()
    for p in db.get("papers", []):
        fh = p.get("file_hash")
        if fh:
            existing_hashes.add(fh)

    pdf_files = sorted(
        [
            os.path.join(pdf_dir, f)
            for f in os.listdir(pdf_dir)
            if f.lower().endswith(".pdf")
        ]
    )

    if not pdf_files:
        print(f"No PDF files found in {pdf_dir}")
        return

    print(f"{'=' * 60}")
    print("PDF DIRECTORY INGESTION")
    print(f"{'=' * 60}")
    print(f"Source dir: {pdf_dir}")
    print(f"PDFs found: {len(pdf_files)}")
    print(f"Database:   {db_path} ({len(db['papers'])} existing papers)")
    print(f"Skip existing: {args.pdf_skip_existing}")
    print()

    new_count = 0
    skipped_count = 0
    for pdf_path in pdf_files:
        fname = os.path.basename(pdf_path)

        file_hash = hash_file(pdf_path)
        if args.pdf_skip_existing and file_hash in existing_hashes:
            print(f"  [SKIP] {fname} (already in database)")
            skipped_count += 1
            continue

        print(f"  [{new_count + 1}] {fname}")

        pdf_text = extract_pdf_text(pdf_path, max_chars=args.pdf_max_chars)
        if not pdf_text:
            print(f"    Could not extract text, skipping")
            continue

        summary, findings, keywords = summarize_paper(
            fname, pdf_path, "", pdf_text, deepseek_key
        )
        if not summary:
            print(f"    Summarization failed, storing raw text excerpt")
            summary = pdf_text[:300]

        entry = {
            "id": file_hash[:12],
            "title": fname,
            "authors": [],
            "abstract": "",
            "source_url": pdf_path,
            "doi": None,
            "arxiv_id": None,
            "pubmed_id": None,
            "pdf_path": pdf_path,
            "file_hash": file_hash,
            "keywords": keywords,
            "summary": summary,
            "key_findings": findings,
            "search_query": "",
            "date_added": datetime.now().strftime("%Y-%m-%d"),
        }

        if add_paper_to_db(db, entry):
            new_count += 1
            existing_hashes.add(file_hash)
            print(f"    Added to database")
        else:
            print(f"    Already in database (duplicate ID)")

    # Record search history entry for this ingestion
    db["search_history"].append(
        {
            "query": f"pdf-dir:{pdf_dir}",
            "date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "results_count": len(pdf_files),
            "new_papers": new_count,
        }
    )

    save_knowledge_db(db, db_path)
    print()
    print(f"{'=' * 60}")
    print(f"Ingestion complete.")
    print(f"  Total PDFs in dir: {len(pdf_files)}")
    print(f"  Newly added:       {new_count}")
    print(f"  Skipped (cached):  {skipped_count}")
    print(f"  Database entries:  {len(db['papers'])}")
    print(f"  Database:          {db_path}")
    print(f"{'=' * 60}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    # Load .env file from repo root before anything else
    load_env_file()
    args = parse_args()

    if args.pdf_dir:
        mode_pdf_dir(args)
    elif args.agent:
        mode_agent(args)
    elif args.query_db:
        mode_query_db(args)
    elif args.query:
        mode_search_and_download(args)
    else:
        print(
            "Error: Must specify one of --query, --query-db, or --agent.",
            file=sys.stderr,
        )
        print("Use --help for usage information.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
