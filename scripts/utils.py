"""
utils.py  --  Shared helpers for the ALSW HiTOM amplicon-seq pipeline.

All functions are pure Python stdlib (no external dependencies).
"""

import csv
import gzip
import json
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# FASTQ I/O
# ---------------------------------------------------------------------------

def open_gz(path, mode="rt"):
    """Open a gzip-compressed or plain-text file."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def fastq_iter(handle):
    """
    Yield (header, seq, plus, qual) tuples from an open FASTQ handle.
    All four fields are stripped of trailing newlines.
    Stops cleanly at EOF; raises ValueError for truncated records.
    """
    while True:
        h = handle.readline()
        if not h:
            return
        s = handle.readline()
        p = handle.readline()
        q = handle.readline()
        if not q:
            raise ValueError("Truncated FASTQ record (incomplete final record)")
        yield (h.rstrip("\n"), s.rstrip("\n"), p.rstrip("\n"), q.rstrip("\n"))


def count_fastq_reads(path):
    """Count the number of read records in a FASTQ.gz file."""
    n = 0
    with open_gz(path) as fh:
        for _ in fastq_iter(fh):
            n += 1
    return n


# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence (case-preserved → uppercase)."""
    return seq.translate(_COMP)[::-1].upper()


def hamming_distance(a, b):
    """Hamming distance between two equal-length strings. Returns inf if lengths differ."""
    if len(a) != len(b):
        return float("inf")
    return sum(x != y for x, y in zip(a, b))


def motif_hit(seq, motif, max_mismatches=0):
    """
    Return True if motif (or its reverse complement) appears in seq
    with at most max_mismatches substitutions.

    For max_mismatches=0 (default, validated pipeline behaviour) this is
    equivalent to a plain substring search and is O(n).
    """
    seq    = seq.upper()
    motif  = motif.upper()
    mot_rc = reverse_complement(motif)
    L      = len(motif)

    if max_mismatches == 0:
        return (motif in seq) or (mot_rc in seq)

    # Mismatch-tolerant search (used only when explicitly requested)
    for target in (motif, mot_rc):
        for i in range(len(seq) - L + 1):
            if hamming_distance(seq[i : i + L], target) <= max_mismatches:
                return True
    return False


def search_motif_in_readpair(s1, s2, motif, max_mismatches=0):
    """
    Search for motif across R1 + R2 + their reverse complements.
    This is the exact logic from the validated run_all_motif_counts.py.

    Combined string: s1 + N + s2 + N + rc(s1) + N + rc(s2)
    """
    combined = (
        s1.upper() + "N"
        + s2.upper() + "N"
        + reverse_complement(s1) + "N"
        + reverse_complement(s2)
    )
    return motif_hit(combined, motif, max_mismatches)


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def load_json(path):
    """Load and return parsed JSON from a file path."""
    with open(path) as fh:
        return json.load(fh)


def load_protocol_constants(path=None):
    """
    Load Hi-TOM protocol constants.
    If path is None, uses pipeline/config/hitom_protocol_constants.json
    relative to this script.
    """
    if path is None:
        path = Path(__file__).resolve().parent.parent / "config" / "hitom_protocol_constants.json"
    return load_json(path)


def load_targets(path=None):
    """
    Load targets.json config.
    If path is None, uses pipeline/config/targets.json relative to this script.
    """
    if path is None:
        path = Path(__file__).resolve().parent.parent / "config" / "targets.json"
    return load_json(path)


# ---------------------------------------------------------------------------
# Barcode resolution
# ---------------------------------------------------------------------------

def resolve_barcode_seq(name_or_seq, fw_dict, rv_dict, side):
    """
    Resolve a barcode name (e.g. 'F1', 'R-A') or a raw 4-bp sequence to
    an uppercase sequence.

    side: 'forward' | 'reverse'
    Raises ValueError with a clear message if the name is not recognised.
    """
    s = name_or_seq.strip().upper()

    if side == "forward":
        if s in fw_dict:
            return fw_dict[s].upper()
        if all(c in "ACGT" for c in s):
            return s
        known = ", ".join(sorted(fw_dict))
        raise ValueError(
            f"Unknown forward barcode '{name_or_seq}'. "
            f"Known names: {known}. "
            f"You may also pass a raw 4-bp sequence (ACGT only)."
        )
    else:
        if s in rv_dict:
            return rv_dict[s].upper()
        if all(c in "ACGT" for c in s):
            return s
        known = ", ".join(sorted(rv_dict))
        raise ValueError(
            f"Unknown reverse barcode '{name_or_seq}'. "
            f"Known names: {known}. "
            f"You may also pass a raw 4-bp sequence (ACGT only)."
        )


def build_barcode_lookup(metadata_rows, fw_dict, rv_dict):
    """
    Build a {(fw_seq, rv_seq): sample_id} mapping from a list of metadata
    dicts, each containing 'sample_id', 'forward_barcode', 'reverse_barcode'.

    Raises ValueError on unknown barcode names or duplicate barcode pairs.
    """
    lookup = {}
    errors = []
    for row in metadata_rows:
        sid = row["sample_id"]
        try:
            fw = resolve_barcode_seq(row["forward_barcode"], fw_dict, rv_dict, "forward")
            rv = resolve_barcode_seq(row["reverse_barcode"], fw_dict, rv_dict, "reverse")
        except ValueError as exc:
            errors.append(f"  sample '{sid}': {exc}")
            continue
        key = (fw, rv)
        if key in lookup:
            errors.append(
                f"  Duplicate barcode pair ({fw}, {rv}) for samples "
                f"'{lookup[key]}' and '{sid}'"
            )
            continue
        lookup[key] = sid

    if errors:
        msg = "Barcode resolution failed:\n" + "\n".join(errors)
        print(msg, file=sys.stderr)
        raise ValueError(msg)

    return lookup


# ---------------------------------------------------------------------------
# Metadata loading
# ---------------------------------------------------------------------------

def load_metadata(path, fw_dict, rv_dict):
    """
    Load sample_metadata.csv, resolve barcodes, and return
    (metadata_rows, barcode_lookup).

    metadata_rows : list of dicts (one per sample)
    barcode_lookup: {(fw_seq, rv_seq): sample_id}
    """
    rows = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            rows.append(row)

    if not rows:
        raise ValueError(f"Metadata file is empty: {path}")

    lookup = build_barcode_lookup(rows, fw_dict, rv_dict)
    return rows, lookup


# ---------------------------------------------------------------------------
# Input FASTQ discovery
# ---------------------------------------------------------------------------

def locate_fastq_pair(directory):
    """
    Find a single R1 / R2 FASTQ.gz pair in a directory.
    Matches *_R1_*.fastq.gz or *_R1.fastq.gz (and equivalent R2).
    Returns (r1_path, r2_path) as Path objects.
    Raises FileNotFoundError or ValueError if not exactly one pair.
    """
    directory = Path(directory)
    r1_cands = sorted(directory.glob("*_R1_*.fastq.gz")) + \
               sorted(directory.glob("*_R1.fastq.gz"))
    r2_cands = sorted(directory.glob("*_R2_*.fastq.gz")) + \
               sorted(directory.glob("*_R2.fastq.gz"))

    # Deduplicate (a file matching both patterns would appear twice)
    r1_cands = list(dict.fromkeys(r1_cands))
    r2_cands = list(dict.fromkeys(r2_cands))

    if not r1_cands:
        raise FileNotFoundError(f"No *_R1*.fastq.gz found in {directory}")
    if not r2_cands:
        raise FileNotFoundError(f"No *_R2*.fastq.gz found in {directory}")
    if len(r1_cands) > 1:
        raise ValueError(f"Expected exactly one R1 file, found: {r1_cands}")
    if len(r2_cands) > 1:
        raise ValueError(f"Expected exactly one R2 file, found: {r2_cands}")

    return r1_cands[0], r2_cands[0]


def locate_sample_fastqs(demux_clean_dir, sample_ids=None):
    """
    Return a list of (sample_id, r1_path, r2_path) tuples from a
    demux_clean directory.  If sample_ids is given, only those samples
    are returned (in that order).
    """
    demux_clean_dir = Path(demux_clean_dir)
    results = []
    candidates = sorted(demux_clean_dir.glob("*_R1.fastq.gz"))

    for r1 in candidates:
        sid = r1.name.replace("_R1.fastq.gz", "")
        if sample_ids is not None and sid not in sample_ids:
            continue
        r2 = demux_clean_dir / f"{sid}_R2.fastq.gz"
        if r2.exists():
            results.append((sid, r1, r2))

    return results


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_csv(rows, path, fieldnames=None):
    """Write a list-of-dicts to a CSV file.  Creates parent dirs as needed."""
    if not rows:
        return
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def write_tsv(rows, path, fieldnames=None):
    """Write a list-of-dicts to a TSV file.  Creates parent dirs as needed."""
    if not rows:
        return
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        # lineterminator="\n" produces Unix LF endings (not CRLF default from RFC 4180)
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def read_tsv(path):
    """Read a TSV file and return a list of dicts."""
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_csv_file(path):
    """Read a CSV file and return a list of dicts."""
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))
