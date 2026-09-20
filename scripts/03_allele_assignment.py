#!/usr/bin/env python3
"""
03_allele_assignment.py  --  Allele-specific read assignment using discriminating motifs.

For targets with validated allele-discriminating flanking SNPs (currently ALSP
chr01 vs chr11), counts reads matching the chromosome-specific motif context.

For ALSW: the validated analysis uses a single WT/desired motif pair with NO
chr01/chr11 allele split.  ALSW reads are therefore reported as 'single_allele'.
DO NOT add ALSW allele discrimination unless validated flanking motifs are known.

Logic ported from: run_all_motif_counts.py (the validated pipeline).
Motif search uses the same combined-strand search:
  combined = R1 + N + R2 + N + RC(R1) + N + RC(R2)

Outputs (all under <output>/03_allele/):
  allele_counts.tsv          per-target per-sample allele counts (WT + desired)
  allele_specific_summary.csv summary focused on chr01/chr11 split
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils

PIPELINE_ROOT    = Path(__file__).resolve().parent.parent
DEFAULT_METADATA = PIPELINE_ROOT / "config" / "sample_metadata.csv"
DEFAULT_CONFIG   = PIPELINE_ROOT / "config" / "targets.json"


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input",    required=True,
                   help="Root output directory (expects <input>/01_demux/demux_clean/)")
    p.add_argument("--output",   required=True,
                   help="Root output directory (outputs go into <output>/03_allele/)")
    p.add_argument("--metadata", default=str(DEFAULT_METADATA))
    p.add_argument("--config",   default=str(DEFAULT_CONFIG),
                   help="targets.json")
    return p.parse_args()


def main():
    args      = parse_args()
    clean_dir = Path(args.input) / "01_demux" / "demux_clean"
    out_dir   = Path(args.output) / "03_allele"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not clean_dir.exists():
        print(f"[03_allele] ERROR: demux_clean not found: {clean_dir}", file=sys.stderr)
        sys.exit(1)

    meta_rows = utils.read_csv_file(args.metadata)
    # Build sample_id → target mapping
    sample_to_target = {r["sample_id"]: r["target"] for r in meta_rows}

    targets_config = utils.load_targets(args.config)
    targets        = targets_config["targets"]

    # Only process targets that have an allele field (chr01/chr11 etc.)
    # Plus single-allele targets (allele = null) for completeness
    # Filter to only targets categorised as allele_and_edit or edit_classification
    allele_targets = [t for t in targets
                      if t.get("category") in ("allele_and_edit", "edit_classification")]

    # Iterate in metadata order (preserves sample ordering from validated analysis)
    samples = []
    for row in meta_rows:
        sid = row["sample_id"]
        r1  = clean_dir / f"{sid}_R1.fastq.gz"
        r2  = clean_dir / f"{sid}_R2.fastq.gz"
        if r1.exists() and r2.exists():
            samples.append((sid, r1, r2))

    counts_rows   = []
    summary_rows  = []

    total_targets = len(allele_targets)
    total_samples = len(samples)
    print(f"[03_allele] {total_targets} target motifs × {total_samples} samples")

    for target in allele_targets:
        tname         = target["name"]
        wt_motif      = target["wt_motif"]
        desired_motif = target["desired_motif"]
        applicable    = [a.upper() for a in target["applicable_to"]]
        allele_label  = target.get("allele") or "single_allele"

        for sid, r1_path, r2_path in samples:
            # Skip samples not applicable to this target
            sample_target = sample_to_target.get(sid, "").upper()
            if sample_target not in applicable:
                continue

            total = wt_hits = desired_hits = both_hits = neither_hits = 0

            with utils.open_gz(r1_path) as f1, utils.open_gz(r2_path) as f2:
                for (_, s1, _, _), (_, s2, _, _) in zip(
                        utils.fastq_iter(f1), utils.fastq_iter(f2)):
                    has_wt      = utils.search_motif_in_readpair(s1, s2, wt_motif)
                    has_desired = utils.search_motif_in_readpair(s1, s2, desired_motif)
                    total += 1
                    if has_wt:      wt_hits     += 1
                    if has_desired: desired_hits += 1
                    if has_wt and has_desired: both_hits    += 1
                    if not has_wt and not has_desired: neither_hits += 1

            informative = wt_hits + desired_hits
            desired_pct = (desired_hits / informative * 100) if informative else 0.0

            counts_rows.append({
                "target":             tname,
                "allele":             allele_label,
                "sample":             sid,
                "total_read_pairs":   total,
                "wt_motif_hits":      wt_hits,
                "desired_motif_hits": desired_hits,
                "both_hits":          both_hits,
                "neither_hits":       neither_hits,
                "informative_reads":  informative,
                "desired_pct_among_motif_hits": f"{desired_pct:.4f}",
            })

            # Collect chr01/chr11 pairs for summary
            if allele_label in ("chr01", "chr11"):
                summary_rows.append({
                    "target_group":   tname.split("_chr")[0] if "_chr" in tname else tname,
                    "allele":         allele_label,
                    "sample":         sid,
                    "wt_hits":        wt_hits,
                    "desired_hits":   desired_hits,
                    "desired_pct":    f"{desired_pct:.4f}",
                })

    utils.write_tsv(counts_rows, out_dir / "allele_counts.tsv")
    utils.write_csv(summary_rows, out_dir / "allele_specific_summary.csv")

    print(f"[03_allele] Wrote {len(counts_rows)} rows to allele_counts.tsv")
    print(f"[03_allele] Done. All outputs in {out_dir}")


if __name__ == "__main__":
    main()
