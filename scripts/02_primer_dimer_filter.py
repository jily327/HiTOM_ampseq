#!/usr/bin/env python3
"""
02_primer_dimer_filter.py  --  Read-length QC and primer-dimer flagging.

*** NEW QC CODE — not part of the validated motif-count analysis. ***

This script reads the post-demux, tail-bridge-clipped FASTQs from 01_demux
and computes read-length distributions.  Short reads (< --min-amplicon-len,
default 100 bp) are flagged as potential primer dimers or short non-specific
products.

A notable feature of this dataset: the R1 read-length distribution shows a
peak around 87 bp which corresponds to short amplification products / primer
dimers (observed during the original analysis session on 2026-06-01).

IMPORTANT: This script is INFORMATIONAL ONLY.
  - It does NOT remove reads from the demux_clean FASTQs.
  - Downstream scripts (03, 04) still use all reads from demux_clean/.
  - The primer-dimer fraction reported here is a QC metric only.
  - If you want to filter short reads before motif counting, you must
    do so explicitly and re-run 03/04 on the filtered FASTQs.

Outputs (all under <output>/02_qc/):
  primer_dimer_qc_summary.csv        per-sample short-read counts + fractions
  per_sample_length_histogram.csv    R1 read-length histograms (1-bp bins)
  potential_dimer_reads_summary.csv  reads below threshold by sample
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils

PIPELINE_ROOT    = Path(__file__).resolve().parent.parent
DEFAULT_METADATA = PIPELINE_ROOT / "config" / "sample_metadata.csv"

MIN_AMPLICON_LEN_DEFAULT = 100   # reads shorter than this are flagged


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input",  required=True,
                   help="Root output directory from run_all.sh "
                        "(expects <input>/01_demux/demux_clean/)")
    p.add_argument("--output", required=True,
                   help="Root output directory (outputs go into <output>/02_qc/)")
    p.add_argument("--metadata", default=str(DEFAULT_METADATA))
    p.add_argument("--min-amplicon-len", type=int, default=MIN_AMPLICON_LEN_DEFAULT,
                   dest="min_amplicon_len",
                   help=f"R1 length threshold for primer-dimer flag "
                        f"[default: {MIN_AMPLICON_LEN_DEFAULT}]")
    return p.parse_args()


def main():
    args     = parse_args()
    clean_dir = Path(args.input) / "01_demux" / "demux_clean"
    out_dir   = Path(args.output) / "02_qc"
    out_dir.mkdir(parents=True, exist_ok=True)
    min_len   = args.min_amplicon_len

    if not clean_dir.exists():
        print(f"[02_qc] ERROR: demux_clean directory not found: {clean_dir}", file=sys.stderr)
        print(f"[02_qc]   Run 01_demux_and_qc.py first.", file=sys.stderr)
        sys.exit(1)

    print(f"[02_qc] Reading from {clean_dir}")
    print(f"[02_qc] Short-read threshold: < {min_len} bp (R1)")

    # Build list of samples from metadata (skip undetermined)
    meta_rows = utils.read_csv_file(args.metadata)
    sample_ids = {r["sample_id"] for r in meta_rows}

    samples = utils.locate_sample_fastqs(clean_dir, sample_ids=None)

    qc_rows   = []
    hist_rows = []
    short_rows = []

    for sid, r1_path, r2_path in samples:
        lengths = []
        with utils.open_gz(r1_path) as fh:
            for _, s, _, _ in utils.fastq_iter(fh):
                lengths.append(len(s))

        if not lengths:
            continue

        total         = len(lengths)
        short_count   = sum(1 for l in lengths if l < min_len)
        long_count    = total - short_count
        short_pct     = short_count / total * 100 if total else 0

        # Collect the 87-bp peak specifically (±5 bp window)
        peak87_count  = sum(1 for l in lengths if 82 <= l <= 92)
        peak87_pct    = peak87_count / total * 100 if total else 0

        qc_rows.append({
            "sample":                      sid,
            "total_read_pairs":            total,
            "short_reads_lt_threshold":    short_count,
            "short_read_threshold_bp":     min_len,
            "short_fraction_pct":          f"{short_pct:.2f}",
            "potential_amplicon_reads":    long_count,
            "amplicon_fraction_pct":       f"{100 - short_pct:.2f}",
            "reads_82_92bp_peak":          peak87_count,
            "reads_82_92bp_pct":           f"{peak87_pct:.2f}",
            "QC_NOTE":                     "INFORMATIONAL_ONLY_reads_not_filtered",
        })

        if short_count > 0:
            short_rows.append({
                "sample":                 sid,
                "short_reads_lt_threshold": short_count,
                "threshold_bp":           min_len,
                "short_fraction_pct":     f"{short_pct:.2f}",
            })

        # Length histogram (bin by bp)
        hist = defaultdict(int)
        for l in lengths:
            hist[l] += 1
        for bp in sorted(hist):
            hist_rows.append({"sample": sid, "R1_length_bp": bp, "count": hist[bp]})

    utils.write_csv(qc_rows,   out_dir / "primer_dimer_qc_summary.csv")
    utils.write_csv(hist_rows, out_dir / "per_sample_length_histogram.csv")
    utils.write_csv(short_rows or [{"note": "No samples exceeded short-read threshold"}],
                    out_dir / "potential_dimer_reads_summary.csv")

    # Summary print
    total_samples = len(qc_rows)
    high_dimer = [r for r in qc_rows if float(r["short_fraction_pct"]) > 20]
    print(f"[02_qc] Processed {total_samples} samples")
    if high_dimer:
        print(f"[02_qc] WARNING: {len(high_dimer)} samples have >20% short reads:")
        for r in high_dimer:
            print(f"  {r['sample']}: {r['short_fraction_pct']}% short")
    print(f"[02_qc] Done. All outputs in {out_dir}")
    print(f"[02_qc] NOTE: Reads were NOT filtered. "
          f"These are QC metrics only; all reads pass to downstream scripts.")


if __name__ == "__main__":
    main()
