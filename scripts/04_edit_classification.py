#!/usr/bin/env python3
"""
04_edit_classification.py  --  WT vs desired edit motif counting.

*** VALIDATED PIPELINE STEP ***

Faithful port of run_all_motif_counts.py.  Produces motif_count_summary.tsv
which is the gold-standard output used for all biological conclusions.

For each target × applicable-sample combination, counts:
  - total read pairs
  - wt_motif_hits
  - desired_motif_hits
  - both_hits
  - neither_hits
  - informative_reads = wt_hits + desired_hits
  - desired_percent_among_motif_hits = desired_hits / informative * 100

Motif search logic (identical to validated run_all_motif_counts.py):
  combined = R1 + N + R2 + N + RC(R1) + N + RC(R2)
  hit = (motif in combined) or (RC(motif) in combined)

All 7 target groups are processed:
  ALSP_core, ALSP_chr01, ALSP_chr11  (ALSP samples only)
  ALSW                                (ALSW samples only)
  EPSPS_full, EPSPS_partial1_front, EPSPS_partial2_back  (EPSPS samples only)

NOT-VALIDATED columns (stubs):
  indel_reads            = NOT_VALIDATED_alignment_required
  imprecise_PE_reads     = NOT_VALIDATED_alignment_required
  Indel and imprecise-PE classification requires alignment-based calling
  (e.g. CRISPResso2 or BWA + variant calling). These columns are zero-filled
  stubs so the table schema is forward-compatible.  See README.md section
  "Known limitations".

Outputs (all under <output>/04_edit/):
  edit_counts.tsv          full table including stub columns
  motif_count_summary.tsv  gold-standard format (matches validated output exactly)
"""

import argparse
import sys
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
                   help="Root output directory (outputs go into <output>/04_edit/)")
    p.add_argument("--metadata", default=str(DEFAULT_METADATA))
    p.add_argument("--config",   default=str(DEFAULT_CONFIG))
    return p.parse_args()


def main():
    args      = parse_args()
    clean_dir = Path(args.input) / "01_demux" / "demux_clean"
    out_dir   = Path(args.output) / "04_edit"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not clean_dir.exists():
        print(f"[04_edit] ERROR: demux_clean not found: {clean_dir}", file=sys.stderr)
        sys.exit(1)

    meta_rows = utils.read_csv_file(args.metadata)
    sample_to_target = {r["sample_id"]: r["target"] for r in meta_rows}

    targets_config = utils.load_targets(args.config)
    targets        = targets_config["targets"]

    # Build sample list in metadata order (preserves row order of gold-standard output)
    samples = []
    for row in meta_rows:
        sid = row["sample_id"]
        r1  = clean_dir / f"{sid}_R1.fastq.gz"
        r2  = clean_dir / f"{sid}_R2.fastq.gz"
        if r1.exists() and r2.exists():
            samples.append((sid, r1, r2))

    total_targets = len(targets)
    total_samples = len(samples)
    print(f"[04_edit] {total_targets} target motifs × {total_samples} samples")

    edit_rows      = []   # full table with stub columns
    motif_rows     = []   # gold-standard motif_count_summary.tsv format

    for target in targets:
        tname         = target["name"]
        wt_motif      = target["wt_motif"]
        desired_motif = target["desired_motif"]
        applicable    = [a.upper() for a in target["applicable_to"]]

        for sid, r1_path, r2_path in samples:
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
                    if has_wt:      wt_hits      += 1
                    if has_desired: desired_hits  += 1
                    if has_wt and has_desired:          both_hits    += 1
                    if not has_wt and not has_desired:  neither_hits += 1

            informative = wt_hits + desired_hits
            desired_pct = (desired_hits / informative * 100) if informative else 0.0

            # ---------------------------------------------------------------
            # motif_count_summary.tsv  (gold-standard exact format)
            # ---------------------------------------------------------------
            motif_rows.append({
                "target":                          tname,
                "sample":                          sid,
                "total_read_pairs":                total,
                "wt_motif_hits":                   wt_hits,
                "desired_motif_hits":              desired_hits,
                "both_hits":                       both_hits,
                "neither_hits":                    neither_hits,
                "informative_reads":               informative,
                "desired_percent_among_motif_hits": f"{desired_pct:.4f}",
            })

            # ---------------------------------------------------------------
            # edit_counts.tsv  (extended table with NOT_VALIDATED stubs)
            # NOT_VALIDATED columns require alignment-based classification.
            # See README.md "Known limitations".
            # ---------------------------------------------------------------
            edit_rows.append({
                "target":                          tname,
                "sample":                          sid,
                "total_read_pairs":                total,
                "wt_motif_hits":                   wt_hits,
                "desired_motif_hits":              desired_hits,
                "both_hits":                       both_hits,
                "neither_hits":                    neither_hits,
                "informative_reads":               informative,
                "desired_percent_among_motif_hits": f"{desired_pct:.4f}",
                # --- NOT_VALIDATED stubs ---
                "indel_reads":                     "NOT_VALIDATED_alignment_required",
                "imprecise_PE_reads":              "NOT_VALIDATED_alignment_required",
            })

    utils.write_tsv(motif_rows, out_dir / "motif_count_summary.tsv")
    utils.write_tsv(edit_rows,  out_dir / "edit_counts.tsv")

    print(f"[04_edit] Wrote {len(motif_rows)} rows to motif_count_summary.tsv")
    print(f"[04_edit] Done. All outputs in {out_dir}")


if __name__ == "__main__":
    main()
