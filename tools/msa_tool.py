#!/usr/bin/env python3
"""
MSA Generation Tool
===================

Generate multiple sequence alignments (MSA) using jackhmmer (via qjackhmmer)
against a user-specified or default sequence database.  For each matched
sequence, compute its percent identity (sequence similarity) to the input
reference sequence.

The tool produces:
  1. A raw A3M file from jackhmmer
  2. An A3M file after filtering (dedup, length-filter)
  3. A JSON report listing every hit with its similarity score and sequences
  4. A plain-text similarity report for quick inspection
  5. A CSV file with hits above a similarity cutoff (default 0.3)
  6. A FASTA file with full-length sequences (Chai-1 format) for hits
     above the cutoff and within ±15 aa of the query length

Tools and default database are defined in data/config.json under "msa".

Usage:
  python tools/msa_tool.py --fasta query.fasta --output ./msa_results
  python tools/msa_tool.py --fasta query.fasta --output ./msa_results \\
    --db /path/to/custom_db.fasta
  python tools/msa_tool.py --fasta query.fasta --output ./msa_results \\
    --incE 1e-5 --n-iter 5
  python tools/msa_tool.py --fasta query.fasta --output ./msa_results \\
    --sim-cutoff 0.5 --len-tolerance 30

References:
  - qjackhmmer is HMMER 3.1b2 jackhmmer wrapper
    http://hmmer.org/
"""

import argparse
import json
import os
import shutil
import subprocess as sp
import sys
import uuid

from _config import load_config


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
MSA_DELETION_AA = "acdefghiklmnpqrstvwyx"
MSA_RESIDUE_AA = "ACDEFGHIKLMNPQRSTVWYX"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run_command(cmd, desc=None):
    """Run a shell command, print status, and exit on failure."""
    tag = desc or cmd
    print(f"[RUN] {tag}")
    ret = sp.run(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                 universal_newlines=True)
    if ret.returncode != 0:
        print(f"[ERROR] Command failed (exit {ret.returncode}): {cmd}",
              file=sys.stderr)
        if ret.stderr:
            print(ret.stderr, file=sys.stderr)
        sys.exit(1)


def read_fasta(fasta_path):
    """Return (header, seq) for the first entry in a FASTA file."""
    seq_lines = []
    header = None
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is None:
                    header = line[1:]  # strip '>'
                continue
            seq_lines.append(line)
    if header is None:
        print(f"[ERROR] No FASTA header found in {fasta_path}", file=sys.stderr)
        sys.exit(1)
    return header, "".join(seq_lines).replace(" ", "").upper()


def parse_a3m(a3m_path, max_sequences=10000):
    """Parse an A3M file and return sequence info.

    Returns
    -------
    list of dict, each with keys:
        seq_id       – sequence name from header
        align_seq    – aligned sequence (deletion chars removed)
        clean_seq    – aligned sequence with gaps removed
        origin_seq   – original (unaligned) sequence
    """
    hits = []
    seq_id = ""
    seq = ""
    oseq = ""
    count = 0

    with open(a3m_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if count >= max_sequences:
                    break
                if seq and seq_id:
                    hits.append({
                        "seq_id": seq_id,
                        "align_seq": seq,
                        "clean_seq": "".join(
                            c for c in seq.replace("X", "-") if c != "-"
                        ),
                        "origin_seq": oseq,
                    })
                seq_id = line.strip("\n")[1:].replace(",", "_")
                seq = ""
                oseq = ""
                count += 1
            else:
                seq += "".join(c for c in line.strip("\n")
                               if c not in MSA_DELETION_AA)
                oseq += line.split()[0]

    # last entry
    if seq and seq_id:
        hits.append({
            "seq_id": seq_id,
            "align_seq": seq,
            "clean_seq": "".join(
                c for c in seq.replace("X", "-") if c != "-"
            ),
            "origin_seq": oseq,
        })

    return hits


def compute_similarities(hits, ref_seq):
    """Compute per-residue identity to the reference sequence for each hit.

    The reference sequence is the aligned (gap-stripped) wild-type sequence.
    Each hit's alignment column is compared to the matching column in ref_seq;
    only columns where both have an amino acid (not gap) are compared.
    Percent identity = matches / (matches + mismatches).

    Modifies each hit dict in-place by adding:
        "identity" – float between 0.0 and 1.0
    """
    for hit in hits:
        aln = hit["align_seq"]
        if len(aln) != len(ref_seq):
            hit["identity"] = 0.0
            continue
        matches = 0
        total = 0
        for a, r in zip(aln, ref_seq):
            if a != "-" and r != "-":
                total += 1
                if a.upper() == r.upper():
                    matches += 1
        hit["identity"] = matches / total if total > 0 else 0.0

    return hits


# ---------------------------------------------------------------------------
# Full-length sequence extraction from DB
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Full-length sequence extraction from DB
# ---------------------------------------------------------------------------
# Threshold: DBs smaller than this size (bytes) are loaded into memory.
_CUSTOM_DB_SIZE_LIMIT = 500 * 1024 * 1024  # 500 MB




def _extract_db_id(seq_id):
    """Extract a searchable DB identifier from an A3M sequence header.

    HMMER A3M headers look like::

        UniRef90_A0A5A9P0L4/1-123  ... RepID=A0A5A9P0L4_9TELE

    or, for custom DBs, simply whatever token was after '>' in the DB.

    Returns the first whitespace-delimited token (minus any ``/range``
    suffix), which is also the token used as the FASTA header key.
    """
    # Strip /range suffix (e.g. "UniRef90_A0A5A9P0L4/1-123")
    idx = seq_id.find("/")
    stem = seq_id[:idx] if idx != -1 else seq_id
    # Split on whitespace and take first token
    stem = stem.split()[0] if " " in stem else stem
    return stem


def _load_fasta_into_dict(db_path):
    """Load an entire small FASTA file into ``{header_token: full_seq}``.

    Used for custom (user-defined) DBs that are small enough to fit in
    memory.  The ``header_token`` is the first whitespace-delimited token
    after ``>``.
    """
    print(f"[INFO] Loading custom DB into memory: {db_path}")
    lookup = {}
    current_key = None
    seq_parts = []
    with open(db_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_key and seq_parts:
                    lookup[current_key] = "".join(seq_parts).upper()
                current_key = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
    if current_key and seq_parts:
        lookup[current_key] = "".join(seq_parts).upper()
    print(f"[INFO] Loaded {len(lookup)} entries from custom DB")
    return lookup


def _build_fasta_index(db_path, index_path=None):
    """Build a byte-offset index for a large FASTA.

    For very large DBs (e.g. UniRef90 > 1 GB), we build a lightweight
    on-disk index at ``{db_path}.msa_idx.json`` that maps
    ``header_token -> {"offset": int, "len": int}``.  The index is
    cached to avoid re-scanning on subsequent runs.
    """
    if index_path is None:
        index_path = db_path + ".msa_idx.json"

    if os.path.exists(index_path):
        with open(index_path) as fh:
            return json.load(fh)

    print(f"[INFO] Building FASTA index for {db_path} (this is done once)")
    idx = {}
    with open(db_path) as fh:
        offset = 0
        line = fh.readline()
        while line:
            if line.startswith(">"):
                key = line[1:].split()[0]
                idx[key] = {"offset": offset}
            offset = fh.tell()
            line = fh.readline()

    # Second pass: compute byte-length for each entry
    keys = list(idx.keys())
    offsets = sorted((idx[k]["offset"], k) for k in keys)
    for i, (off, key) in enumerate(offsets):
        next_off = offsets[i + 1][0] if i + 1 < len(offsets) else \
            os.path.getsize(db_path)
        idx[key]["len"] = next_off - off

    with open(index_path, "w") as fh:
        json.dump(idx, fh)

    print(f"[INFO] Index written to {index_path} ({len(idx)} entries)")
    return idx


def _extract_full_from_index(seq_id, db_path, db_index):
    """Seek into the large DB using the byte-offset index.

    Returns ``(full_header_str, full_seq_str)`` or ``(None, None)``.
    """
    token = _extract_db_id(seq_id)
    entry = db_index.get(token)
    if entry is None:
        return None, None

    seq_header = None
    seq_parts = []
    end_offset = entry["offset"] + entry["len"]
    with open(db_path) as fh:
        fh.seek(entry["offset"])
        line = fh.readline()
        while line:
            if line.startswith(">"):
                if seq_header is None:
                    seq_header = line[1:].strip()
            else:
                if line.strip():
                    seq_parts.append(line.strip())
            # Use file.tell() after readline() — safe because
            # readline() returns position, while implicit next()
            # from "for line in fh" disables tell().
            if fh.tell() >= end_offset:
                break
            line = fh.readline()

    full_seq = "".join(seq_parts).upper()
    return seq_header, full_seq


def get_full_sequence_lookup(db_path):
    """Return a callable that extracts full-length sequences from *db_path*.

    Inspects the DB file size to decide the strategy:

    * **Small DB** (<= 500 MB) — load entire file into an in-memory dict.
      Ideal for user-defined custom DBs.
    * **Large DB** (> 500 MB) — build a lightweight byte-offset index
      (cached to disk) and seek on demand.  Suitable for UniRef90 etc.

    The returned function has signature ``(seq_id: str) -> (header, seq)``.
    """
    size = os.path.getsize(db_path)
    if size <= _CUSTOM_DB_SIZE_LIMIT:
        lookup = _load_fasta_into_dict(db_path)

        def _lookup(seq_id):
            token = _extract_db_id(seq_id)
            full_seq = lookup.get(token)
            if full_seq is None:
                return None, None
            return token, full_seq
        return _lookup
    else:
        db_index = _build_fasta_index(db_path)

        def _lookup(seq_id):
            return _extract_full_from_index(seq_id, db_path, db_index)
        return _lookup


def _build_seq_to_header_map(raw_a3m_path, exclude_header=None):
    """Parse the raw A3M and return ``{clean_seq: db_header_token}``.

    Only keeps entries whose header does NOT match *exclude_header*
    (typically the query header), so subsequent lookups return the
    original DB header rather than the synthetic query header.
    """
    mapping = {}
    if not os.path.exists(raw_a3m_path):
        return mapping
    seq_id = ""
    seq = ""
    with open(raw_a3m_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq and seq_id:
                    clean = "".join(
                        c for c in seq.replace("X", "-")
                        if c not in MSA_DELETION_AA
                    )
                    clean = "".join(c for c in clean if c != "-")
                    if clean not in mapping:
                        mapping[clean] = _extract_db_id(seq_id)
                seq_id = line.strip("\n")[1:].replace(",", "_")
                seq = ""
            else:
                seq += "".join(c for c in line.strip("\n")
                               if c not in MSA_DELETION_AA)
    if seq and seq_id:
        clean = "".join(
            c for c in seq.replace("X", "-") if c not in MSA_DELETION_AA
        )
        clean = "".join(c for c in clean if c != "-")
        if clean not in mapping:
            mapping[clean] = _extract_db_id(seq_id)
    # Remove entries that came from the query itself
    if exclude_header:
        exclude_token = _extract_db_id(exclude_header)
        to_delete = [k for k, v in mapping.items() if v == exclude_token]
        for k in to_delete:
            del mapping[k]
    return mapping


# ---------------------------------------------------------------------------
# Core function
# ---------------------------------------------------------------------------

def generate_msa(fasta_path, output_dir, db_path, inc_e, n_iter, max_seqs,
                 skip_existing, config, sim_cutoff=0.3, len_tolerance=15):
    """Run jackhmmer MSA search, filter, compute similarities, and produce
    CSV + Chai-1 FASTA outputs for high-similarity hits.

    Parameters
    ----------
    fasta_path    : str – input query FASTA
    output_dir    : str – output directory
    db_path       : str – sequence database to search against
    inc_e         : float – inclusion E-value threshold
    n_iter        : int – number of jackhmmer iterations
    max_seqs      : int – maximum sequences to keep in parsed result
    skip_existing : bool – skip if output files exist
    config        : dict – tool config (paths to binaries)
    sim_cutoff    : float – similarity cutoff (0-1) for CSV/FASTA output
    len_tolerance : int – max length difference in aa for FASTA output

    Returns
    -------
    dict with keys:
        raw_a3m, filtered_a3m, report_json, report_txt, csv_path,
        fasta_path, num_hits, num_passed, avg_identity, reference_seq
    """
    os.makedirs(output_dir, exist_ok=True)

    query_name = os.path.splitext(os.path.basename(fasta_path))[0]
    raw_a3m = os.path.join(output_dir, f"{query_name}.a3m")
    filtered_a3m = os.path.join(output_dir, f"{query_name}.filtered.a3m")
    report_json = os.path.join(output_dir, f"{query_name}.msa_report.json")
    report_txt = os.path.join(output_dir, f"{query_name}.similarity.txt")
    csv_path = os.path.join(output_dir, f"{query_name}.msa_hits.csv")
    fasta_path_out = os.path.join(output_dir, f"{query_name}.msa_hits.fasta")

    # -- read reference sequence
    ref_header, ref_seq = read_fasta(fasta_path)
    print(f"[INFO] Query: {ref_header} ({len(ref_seq)} residues)")

    # -- Step 1: run jackhmmer
    if skip_existing and os.path.exists(raw_a3m):
        print(f"[SKIP] {raw_a3m} exists, skipping jackhmmer search")
    else:
        tmp_raw = f"/tmp/msa_{uuid.uuid4().hex}.a3m"
        print(f"[INFO] Running jackhmmer search against {db_path}")
        cmd = (
            f"{config['msa']['qjackhmmer']} "
            f"-N {n_iter} "
            f"-B {tmp_raw} "
            f"--incE {inc_e} "
            f"{fasta_path} {db_path} > /dev/null 2>&1"
        )
        run_command(cmd, "jackhmmer search")
        if os.path.exists(tmp_raw):
            shutil.copy2(tmp_raw, raw_a3m)
        else:
            print(f"[ERROR] jackhmmer failed to produce output at {tmp_raw}",
                  file=sys.stderr)
            sys.exit(1)

    # -- Step 2: filter
    if skip_existing and os.path.exists(filtered_a3m):
        print(f"[SKIP] {filtered_a3m} exists, skipping filtering")
    else:
        print(f"[INFO] Filtering MSA (keeping sequences >= 80% query length)")
        _filter_msa(raw_a3m, filtered_a3m, min_len_ratio=0.8)

    # -- Step 3: parse and compute similarities
    ref_align_seq = _make_reference_align_seq(ref_seq, filtered_a3m)

    hits = parse_a3m(filtered_a3m, max_sequences=max_seqs)
    print(f"[INFO] Parsed {len(hits)} sequences from filtered MSA")
    if not hits:
        print("[WARN] No hits found in MSA output")
        hits = []

    # Ensure query is first hit
    has_query = any(h["seq_id"] == ref_header for h in hits)
    if not has_query:
        hits.insert(0, {
            "seq_id": ref_header,
            "align_seq": ref_align_seq,
            "clean_seq": ref_seq,
            "origin_seq": ref_seq,
        })

    compute_similarities(hits, ref_align_seq)

    # -- Step 4: extract full-length sequences from DB for high-similarity hits
    print(f"[INFO] Building sequence lookup from {db_path}...")
    db_lookup = get_full_sequence_lookup(db_path)
    # Build a reverse map: clean_seq -> db_header from the raw A3M
    # (before dedup).  This is used when the query is synthetic and its
    # header doesn't appear in the DB.
    seq_to_db_header = _build_seq_to_header_map(
        raw_a3m, exclude_header=ref_header)
    seq_id_counts = {}  # for generating unique seq_id in Chai-1 format

    for hit in hits:
        if hit["identity"] >= sim_cutoff:
            db_header, full_seq = db_lookup(hit["seq_id"])
            # Fallback: if header lookup failed and the hit is identical to
            # the query, look up by sequence to find the original DB header.
            if full_seq is None and hit["identity"] >= 0.999:
                orig_header = seq_to_db_header.get(hit["clean_seq"])
                if orig_header:
                    db_header, full_seq = db_lookup(orig_header)
            hit["full_header"] = db_header or hit["seq_id"]
            hit["full_seq"] = full_seq
            if hit["full_seq"]:
                seq_id_counts[hit["full_seq"]] = \
                    seq_id_counts.get(hit["full_seq"], 0) + 1
        else:
            hit["full_header"] = None
            hit["full_seq"] = None

    # -- Step 5: write JSON report (same as before, now with full_seq)
    identities = [h["identity"] for h in hits]
    avg_identity = sum(identities) / len(identities) if identities else 0.0

    summary = {
        "reference": ref_header,
        "reference_length": len(ref_seq),
        "database": db_path,
        "num_hits": len(hits),
        "avg_identity": round(avg_identity, 4),
        "sim_cutoff": sim_cutoff,
        "hits": [
            {
                "seq_id": h["seq_id"],
                "identity": round(h["identity"], 4),
                "clean_seq": h["clean_seq"],
                "full_seq": h.get("full_seq"),
            }
            for h in hits
        ],
    }
    with open(report_json, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[OUT] JSON report → {report_json}")

    # -- Step 6: write similarity TXT
    with open(report_txt, "w") as fh:
        fh.write(f"Reference: {ref_header}\n")
        fh.write(f"Database: {db_path}\n")
        fh.write(f"Total hits: {len(hits)}\n")
        fh.write(f"Average identity: {avg_identity:.4f}\n")
        fh.write("-" * 72 + "\n")
        for h in sorted(hits, key=lambda x: x["identity"], reverse=True):
            fh.write(f"{h['identity']:.4f}  {h['seq_id']}\n")
    print(f"[OUT] Similarity report → {report_txt}")

    # -- Step 7: write CSV (all hits above cutoff)
    csv_rows = [h for h in hits if h["identity"] >= sim_cutoff]
    with open(csv_path, "w") as fh:
        fh.write("seq_id,identity,full_seq_length,full_seq\n")
        for h in sorted(csv_rows, key=lambda x: x["identity"], reverse=True):
            full_seq = h.get("full_seq") or ""
            fh.write(f"{h['seq_id']},{h['identity']:.4f},{len(full_seq)},"
                     f"{full_seq}\n")
    print(f"[OUT] CSV hits (sim>={sim_cutoff}) → {csv_path} ({len(csv_rows)} entries)")

    # -- Step 8: write Chai-1 format FASTA (hits above cutoff AND length OK)
    num_fasta = 0
    ref_len = len(ref_seq)
    written_seqs = set()
    with open(fasta_path_out, "w") as fh:
        # Always include the reference as the first entry
        fh.write(f">protein|{query_name} {ref_header}\n{ref_seq}\n")
        written_seqs.add(ref_seq)
        num_fasta += 1

        for h in sorted(csv_rows, key=lambda x: x["identity"], reverse=True):
            full_seq = h.get("full_seq")
            if not full_seq:
                continue
            # Skip duplicates (e.g. the query identity=1.0 hit found in DB)
            if full_seq in written_seqs:
                continue
            # Length filter: within tolerance of reference
            if abs(len(full_seq) - ref_len) > len_tolerance:
                continue
            # Build unique name
            db_id = _extract_db_id(h["seq_id"])
            seq_name = f"{db_id}_{seq_id_counts.get(full_seq, 1)}"
            fh.write(f">protein|{seq_name} {h['seq_id']}\n{full_seq}\n")
            written_seqs.add(full_seq)
            num_fasta += 1

    print(f"[OUT] Chai-1 FASTA (sim>={sim_cutoff}, len ±{len_tolerance}) "
          f"→ {fasta_path_out} ({num_fasta} entries)")

    return {
        "raw_a3m": raw_a3m,
        "filtered_a3m": filtered_a3m,
        "report_json": report_json,
        "report_txt": report_txt,
        "csv_path": csv_path,
        "fasta_path": fasta_path_out,
        "num_hits": len(hits),
        "num_passed": len(csv_rows),
        "avg_identity": avg_identity,
        "reference_seq": ref_seq,
    }


def _make_reference_align_seq(ref_seq, a3m_path):
    """Build a dummy aligned reference sequence (no gaps, all residues)."""
    return ref_seq.upper()


def _filter_msa(in_a3m, out_a3m, min_len_ratio=0.8):
    """Filter an A3M file: deduplicate and remove short sequences."""
    entries = []
    with open(in_a3m) as fh:
        lines = fh.readlines()

    seq_starts = [i for i, l in enumerate(lines) if l.startswith(">")]
    if not seq_starts:
        with open(out_a3m, "w"):
            pass
        return

    for idx, start in enumerate(seq_starts):
        end = seq_starts[idx + 1] if idx + 1 < len(seq_starts) else len(lines)
        header = lines[start].strip()
        raw_seq = "".join(l.strip() for l in lines[start + 1:end])
        clean = "".join(c for c in raw_seq.replace("X", "-")
                        if c not in MSA_DELETION_AA)
        entries.append((header, clean, raw_seq))

    if not entries:
        with open(out_a3m, "w"):
            pass
        return

    query_len = len(entries[0][1])
    min_len = int(query_len * min_len_ratio)

    seen = set()
    kept = []
    for header, seq, raw_seq in entries:
        dedup_key = "".join(c for c in seq if c != "-")
        if len(dedup_key) < min_len:
            continue
        if dedup_key in seen:
            continue
        seen.add(dedup_key)
        kept.append((header, seq, raw_seq))

    with open(out_a3m, "w") as fh:
        for header, seq, _ in kept:
            fh.write(f"{header}\n{seq}\n")

    print(f"[INFO] Filtered MSA: {len(entries)} → {len(kept)} sequences")
    return kept


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate MSA with jackhmmer and compute sequence similarities",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--fasta", required=True,
                        help="Input query FASTA file")
    parser.add_argument("--output", "-o", required=True,
                        help="Output directory for MSA files and reports")
    parser.add_argument("--db",
                        help="Sequence database to search (overrides "
                             "default_db in config)")
    parser.add_argument("--incE", type=float, default=None,
                        help="Inclusion E-value threshold "
                             "(default: from config, usually 0.001)")
    parser.add_argument("--n-iter", type=int, default=None,
                        help="Number of jackhmmer iterations "
                             "(default: from config, usually 3)")
    parser.add_argument("--max-seqs", type=int, default=10000,
                        help="Maximum sequences to keep in parsed results "
                             "(default: 10000)")
    parser.add_argument("--sim-cutoff", type=float, default=0.3,
                        help="Similarity cutoff (0-1) for CSV/FASTA output "
                             "(default: 0.3)")
    parser.add_argument("--len-tolerance", type=int, default=15,
                        help="Max length difference in aa for FASTA output "
                             "(default: 15)")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip steps whose output files already exist")
    parser.add_argument("--config",
                        help="Path to tool config JSON (default: "
                             "data/config.json in repo root)")
    return parser.parse_args(argv)


def main():
    args = parse_args()
    config = load_config(args.config)

    db_path = args.db or config.get("msa", {}).get("default_db")
    if not db_path:
        print("[ERROR] No database specified. Provide --db or set "
              "msa.default_db in config.json", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(db_path):
        print(f"[ERROR] Database not found: {db_path}", file=sys.stderr)
        sys.exit(1)

    inc_e = args.incE if args.incE is not None else \
        config.get("msa", {}).get("jackhmmer_incE", 0.001)
    n_iter = args.n_iter if args.n_iter is not None else \
        config.get("msa", {}).get("jackhmmer_n_iterations", 3)

    result = generate_msa(
        fasta_path=args.fasta,
        output_dir=args.output,
        db_path=db_path,
        inc_e=inc_e,
        n_iter=n_iter,
        max_seqs=args.max_seqs,
        skip_existing=args.skip_existing,
        config=config,
        sim_cutoff=args.sim_cutoff,
        len_tolerance=args.len_tolerance,
    )

    print(f"\n[DONE] MSA generation complete.")
    print(f"       Hits: {result['num_hits']}")
    print(f"       Passed sim cutoff: {result['num_passed']}")
    print(f"       Avg identity: {result['avg_identity']:.4f}")
    print(f"       Filtered A3M: {result['filtered_a3m']}")
    print(f"       CSV:           {result['csv_path']}")
    print(f"       Chai-1 FASTA:  {result['fasta_path']}")


if __name__ == "__main__":
    main()
