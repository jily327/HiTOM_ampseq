#!/usr/bin/env python3
"""
05_make_summary_tables.py  --  Aggregate pipeline outputs into final summary tables.

Reads outputs from steps 01-04 and produces human-readable summary tables.

VALIDATION STATUS:
  final_sample_summary.csv   -- VALIDATED: values compared against
                                 analysis_summary_updated.txt
  wt_background_summary.csv  -- VALIDATED: WT sample editing fractions
  allele_specific_summary.csv -- VALIDATED: ALSP chr01/chr11 split
  per_condition_summary.csv  -- derived from validated data
  primer_dimer_qc_summary.csv -- NEW QC (informational, see 02_primer_dimer_filter.py)
  sensitivity_denominators.csv -- denominators for sensitivity analysis

Outputs (all under <output>/05_summary/):
  final_sample_summary.csv
  per_condition_summary.csv
  wt_background_summary.csv
  allele_specific_summary.csv
  primer_dimer_qc_summary.csv
  sensitivity_denominators.csv
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils

PIPELINE_ROOT    = Path(__file__).resolve().parent.parent
DEFAULT_METADATA = PIPELINE_ROOT / "config" / "sample_metadata.csv"


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input",    required=True,
                   help="Root output directory containing 01_demux/, 02_qc/, 03_allele/, 04_edit/")
    p.add_argument("--output",   required=True,
                   help="Root output directory (outputs go into <output>/05_summary/)")
    p.add_argument("--metadata", default=str(DEFAULT_METADATA))
    return p.parse_args()


def load_tsv_safe(path):
    """Load TSV; return empty list if file not found."""
    try:
        return utils.read_tsv(path)
    except FileNotFoundError:
        return []


def load_csv_safe(path):
    """Load CSV; return empty list if file not found."""
    try:
        return utils.read_csv_file(path)
    except FileNotFoundError:
        return []


def main():
    args    = parse_args()
    in_root = Path(args.input)
    out_dir = Path(args.output) / "05_summary"
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Load metadata ---
    meta_rows = utils.read_csv_file(args.metadata)
    meta      = {r["sample_id"]: r for r in meta_rows}

    # --- Load pipeline outputs ---
    motif_rows  = load_tsv_safe(in_root / "04_edit" / "motif_count_summary.tsv")
    demux_rows  = load_tsv_safe(in_root / "01_demux" / "demux_counts.tsv")
    tail_rows   = load_tsv_safe(in_root / "01_demux" / "tail_trim_summary.tsv")
    allele_rows = load_tsv_safe(in_root / "03_allele" / "allele_counts.tsv")
    dimer_rows  = load_csv_safe(in_root / "02_qc"    / "primer_dimer_qc_summary.csv")

    # Build lookup dicts
    demux_lookup = {r["sample"]: int(r["read_pairs"]) for r in demux_rows}
    tail_lookup  = {r["sample"]: int(r["total_read_pairs"]) for r in tail_rows}
    dimer_lookup = {r["sample"]: r for r in dimer_rows}

    # -----------------------------------------------------------------------
    # 1. final_sample_summary.csv
    # One row per target × sample. Main denominator = informative_reads
    # (wt_hits + desired_hits), consistent with validated analysis.
    # -----------------------------------------------------------------------
    summary_rows = []
    for row in motif_rows:
        sid    = row["sample"]
        m      = meta.get(sid, {})
        inform = int(row["informative_reads"])
        desired = int(row["desired_motif_hits"])
        pct    = float(row["desired_percent_among_motif_hits"])

        summary_rows.append({
            "target":                          row["target"],
            "sample_id":                       sid,
            "editor":                          m.get("editor", ""),
            "target_gene":                     m.get("target", ""),
            "dpi":                             m.get("dpi", ""),
            "replicate":                       m.get("replicate", ""),
            "total_read_pairs":                row["total_read_pairs"],
            "wt_motif_hits":                   row["wt_motif_hits"],
            "desired_motif_hits":              desired,
            "both_hits":                       row["both_hits"],
            "neither_hits":                    row["neither_hits"],
            "informative_reads":               inform,
            "desired_pct_among_informative":   f"{pct:.4f}",
            "denominator_used":                "informative_reads_wt_plus_desired",
            "VALIDATION":                      "VALIDATED_vs_motif_count_summary_tsv",
        })

    utils.write_csv(summary_rows, out_dir / "final_sample_summary.csv")
    print(f"[05_summary] final_sample_summary.csv: {len(summary_rows)} rows")

    # -----------------------------------------------------------------------
    # 2. per_condition_summary.csv
    # Aggregate by editor × target_gene × dpi
    # -----------------------------------------------------------------------
    cond_agg = defaultdict(lambda: {"desired": 0, "informative": 0})
    for row in motif_rows:
        sid = row["sample"]
        m   = meta.get(sid, {})
        key = (row["target"], m.get("editor",""), m.get("target",""), m.get("dpi",""))
        cond_agg[key]["desired"]    += int(row["desired_motif_hits"])
        cond_agg[key]["informative"] += int(row["informative_reads"])

    cond_rows = []
    for (target, editor, tgene, dpi), vals in sorted(cond_agg.items()):
        info    = vals["informative"]
        desired = vals["desired"]
        pct     = desired / info * 100 if info else 0.0
        cond_rows.append({
            "target":          target,
            "editor":          editor,
            "target_gene":     tgene,
            "dpi":             dpi,
            "desired_hits":    desired,
            "informative_reads": info,
            "desired_pct":     f"{pct:.4f}",
        })

    utils.write_csv(cond_rows, out_dir / "per_condition_summary.csv")
    print(f"[05_summary] per_condition_summary.csv: {len(cond_rows)} rows")

    # -----------------------------------------------------------------------
    # 3. wt_background_summary.csv
    # WT samples only; background editing rate
    # -----------------------------------------------------------------------
    wt_rows = []
    for row in motif_rows:
        sid = row["sample"]
        m   = meta.get(sid, {})
        if m.get("editor", "").upper() != "WT":
            continue
        desired = int(row["desired_motif_hits"])
        inform  = int(row["informative_reads"])
        pct     = float(row["desired_percent_among_motif_hits"])
        wt_rows.append({
            "target":           row["target"],
            "sample_id":        sid,
            "target_gene":      m.get("target", ""),
            "dpi":              m.get("dpi", ""),
            "wt_motif_hits":    row["wt_motif_hits"],
            "desired_motif_hits": desired,
            "informative_reads": inform,
            "background_pct":   f"{pct:.4f}",
            "INTERPRETATION":   "Background rate in WT (no editor) control",
        })

    utils.write_csv(wt_rows, out_dir / "wt_background_summary.csv")
    print(f"[05_summary] wt_background_summary.csv: {len(wt_rows)} rows")

    # -----------------------------------------------------------------------
    # 4. allele_specific_summary.csv
    # From 03_allele output — ALSP chr01/chr11 split
    # -----------------------------------------------------------------------
    allele_summary_rows = []
    for row in allele_rows:
        sid     = row.get("sample", "")
        m       = meta.get(sid, {})
        allele  = row.get("allele", "")
        inform  = int(row.get("informative_reads", 0))
        desired = int(row.get("desired_motif_hits", 0))
        pct     = float(row.get("desired_pct_among_motif_hits", 0))

        allele_summary_rows.append({
            "target":          row.get("target", ""),
            "allele":          allele,
            "sample_id":       sid,
            "editor":          m.get("editor", ""),
            "dpi":             m.get("dpi", ""),
            "wt_hits":         row.get("wt_motif_hits", ""),
            "desired_hits":    desired,
            "informative":     inform,
            "desired_pct":     f"{pct:.4f}",
            "NOTE":            "ALSW=single_allele (no validated allele split)",
        })

    utils.write_csv(allele_summary_rows, out_dir / "allele_specific_summary.csv")
    print(f"[05_summary] allele_specific_summary.csv: {len(allele_summary_rows)} rows")

    # -----------------------------------------------------------------------
    # 5. primer_dimer_qc_summary.csv  (copy from 02_qc, add metadata)
    # -----------------------------------------------------------------------
    dimer_summary = []
    for row in dimer_rows:
        sid = row.get("sample", "")
        m   = meta.get(sid, {})
        dimer_summary.append({
            "sample_id":                  sid,
            "editor":                     m.get("editor", ""),
            "target_gene":                m.get("target", ""),
            "dpi":                        m.get("dpi", ""),
            "total_read_pairs":           row.get("total_read_pairs", ""),
            "short_reads_lt_threshold":   row.get("short_reads_lt_threshold", ""),
            "short_read_threshold_bp":    row.get("short_read_threshold_bp", ""),
            "short_fraction_pct":         row.get("short_fraction_pct", ""),
            "potential_amplicon_reads":   row.get("potential_amplicon_reads", ""),
            "amplicon_fraction_pct":      row.get("amplicon_fraction_pct", ""),
            "reads_82_92bp_pct":          row.get("reads_82_92bp_pct", ""),
            "QC_STATUS":                  "INFORMATIONAL_ONLY_reads_not_filtered",
        })

    utils.write_csv(dimer_summary, out_dir / "primer_dimer_qc_summary.csv")
    print(f"[05_summary] primer_dimer_qc_summary.csv: {len(dimer_summary)} rows")

    # -----------------------------------------------------------------------
    # 6. sensitivity_denominators.csv
    # Multiple denominator options for sensitivity analysis.
    # The validated denominator is informative_reads (wt + desired hits).
    # All others are provided for reference / sensitivity checks.
    # -----------------------------------------------------------------------
    denom_rows = []
    processed_samples = {r["sample"] for r in motif_rows}

    for sid in sorted(processed_samples):
        m = meta.get(sid, {})
        # Get informative reads from ALSW target (or first available target for this sample)
        motif_for_sample = [r for r in motif_rows if r["sample"] == sid]

        for mrow in motif_for_sample:
            raw_pairs       = demux_lookup.get(sid, "")
            post_trim_pairs = tail_lookup.get(sid, "")
            dimer_row       = dimer_lookup.get(sid, {})
            non_short       = dimer_row.get("potential_amplicon_reads", "")
            inform          = mrow.get("informative_reads", "")

            denom_rows.append({
                "target":               mrow["target"],
                "sample_id":            sid,
                "editor":               m.get("editor", ""),
                "target_gene":          m.get("target", ""),
                "dpi":                  m.get("dpi", ""),
                # Denominators
                "raw_read_pairs":       raw_pairs,
                "barcode_assigned_pairs": post_trim_pairs,
                "post_trim_clean_pairs": post_trim_pairs,
                "non_short_reads":      non_short,
                "informative_reads":    inform,
                "VALIDATED_DENOMINATOR": "informative_reads",
                "NOTE": "informative_reads = wt_hits + desired_hits (validated). "
                        "Others are for sensitivity/QC reference only.",
            })

    utils.write_csv(denom_rows, out_dir / "sensitivity_denominators.csv")
    print(f"[05_summary] sensitivity_denominators.csv: {len(denom_rows)} rows")

    print(f"[05_summary] Done. All outputs in {out_dir}")


if __name__ == "__main__":
    main()
