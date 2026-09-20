#!/usr/bin/env python3
"""
10_compare_motif_vs_alignment.py  --  Compare motif-count vs alignment results.

*** OPTIONAL MODULE — NOT PART OF THE VALIDATED MOTIF-COUNT PIPELINE ***
*** ALIGNMENT RESULTS ARE NOT VALIDATED; INTERPRET WITH CAUTION       ***

Joins the validated motif_count_summary.tsv with the alignment-based
alignment_sample_summary.csv (from 07_alignment_edit_classifier.py) and
compares desired edit percentages by sample.

Large discrepancies between motif-based and alignment-based estimates can
indicate:
  - Reads that contain the desired motif but also carry indels elsewhere
    (captured separately by alignment as "desired_plus_indel")
  - Reads with imprecise prime edits near the motif window
  - Method-specific differences (motif window vs read-level classification)

The validated motif-count result is the reference; alignment results are
the comparator.

Outputs (under <output>/10_comparison/):
  motif_vs_alignment_comparison.csv
  discrepancy_report.md
  motif_vs_alignment_scatter.{pdf,png,svg}  (if matplotlib available)
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils

PIPELINE_ROOT    = Path(__file__).resolve().parent.parent
DEFAULT_METADATA = PIPELINE_ROOT / "config" / "sample_metadata.csv"

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input",     required=True,
                   help="Root output dir (expects 04_edit/ and 07_alignment/)")
    p.add_argument("--output",    required=True,
                   help="Root output dir (results in <output>/10_comparison/)")
    p.add_argument("--metadata",  default=str(DEFAULT_METADATA))
    p.add_argument("--threshold", type=float, default=1.0,
                   help="Flag discrepancies > this many pct points [default: 1.0]")
    p.add_argument("--format",    default="pdf,png,svg")
    return p.parse_args()


def main():
    args      = parse_args()
    in_root   = Path(args.input)
    out_dir   = Path(args.output) / "10_comparison"
    out_dir.mkdir(parents=True, exist_ok=True)
    formats   = [f.strip() for f in args.format.split(",")]

    motif_path = in_root / "04_edit" / "motif_count_summary.tsv"
    align_path = in_root / "07_alignment" / "alignment_sample_summary.csv"

    if not motif_path.exists():
        print(f"[10_compare] ERROR: {motif_path} not found. Run 04 first.",
              file=sys.stderr)
        sys.exit(1)
    if not align_path.exists():
        print(f"[10_compare] ERROR: {align_path} not found. Run 07 first.",
              file=sys.stderr)
        sys.exit(1)

    motif_rows = utils.read_tsv(motif_path)
    align_rows = utils.read_csv_file(align_path)
    meta_rows  = utils.read_csv_file(args.metadata)
    meta_lut   = {r["sample_id"]: r for r in meta_rows}

    # Index alignment rows
    align_lut = {}
    for r in align_rows:
        align_lut[(r["target"], r["sample_id"])] = r

    compare_rows = []
    discrepancies = []

    for mr in motif_rows:
        tgt = mr["target"]
        sid = mr["sample"]
        key = (tgt, sid)
        ar  = align_lut.get(key)

        motif_pct = float(mr["desired_percent_among_motif_hits"])

        if ar is None:
            compare_rows.append({
                "target": tgt, "sample_id": sid,
                "motif_desired_pct": f"{motif_pct:.4f}",
                "alignment_desired_pct": "N/A",
                "difference_pct": "N/A",
                "flagged": False,
                "alignment_wt": "",
                "alignment_precise_desired": "",
                "alignment_indel": "",
                "alignment_imprecise_PE": "",
                "VALIDATION_NOTE": "alignment_not_run_for_this_target",
            })
            continue

        align_pct = float(ar["alignment_precise_desired_pct"])
        diff      = abs(motif_pct - align_pct)
        flagged   = diff > args.threshold

        row = {
            "target":                  tgt,
            "sample_id":               sid,
            "editor":                  meta_lut.get(sid, {}).get("editor", ""),
            "dpi":                     meta_lut.get(sid, {}).get("dpi", ""),
            "motif_desired_pct":       f"{motif_pct:.4f}",
            "alignment_desired_pct":   f"{align_pct:.4f}",
            "difference_pct":          f"{diff:.4f}",
            "flagged":                 flagged,
            "alignment_wt":            ar.get("wt", ""),
            "alignment_precise_desired": ar.get("precise_desired", ""),
            "alignment_indel":         ar.get("indel", ""),
            "alignment_imprecise_PE":  ar.get("imprecise_PE", ""),
            "alignment_total_processed": ar.get("total_reads_processed", ""),
            "VALIDATION_NOTE":         "alignment_column_NOT_VALIDATED",
        }
        compare_rows.append(row)
        if flagged:
            discrepancies.append(row)

    # Ensure all rows have the same columns before writing
    ALL_COLS = ["target", "sample_id", "editor", "dpi",
                "motif_desired_pct", "alignment_desired_pct",
                "difference_pct", "flagged",
                "alignment_wt", "alignment_precise_desired",
                "alignment_indel", "alignment_imprecise_PE",
                "alignment_total_processed", "VALIDATION_NOTE"]
    for r in compare_rows:
        for col in ALL_COLS:
            r.setdefault(col, "")

    utils.write_csv(compare_rows, out_dir / "motif_vs_alignment_comparison.csv",
                    fieldnames=ALL_COLS)
    print(f"[10_compare] Wrote {len(compare_rows)} comparison rows")
    print(f"[10_compare] Flagged discrepancies (>{args.threshold}%): "
          f"{len(discrepancies)}")

    # Discrepancy report
    report_lines = [
        "# Motif-count vs Alignment Comparison Report",
        "",
        f"Discrepancy threshold: {args.threshold} percentage points",
        f"Total comparisons: {len(compare_rows)}",
        f"Flagged (|motif − alignment| > threshold): {len(discrepancies)}",
        "",
        "**Note**: Motif-count result is the validated gold standard for this project.",
        "Alignment results are exploratory and NOT validated.",
        "",
    ]
    if discrepancies:
        report_lines += [
            "## Flagged samples",
            "",
            "| Target | Sample | Motif % | Alignment % | Diff |",
            "|--------|--------|---------|-------------|------|",
        ]
        for r in discrepancies:
            report_lines.append(
                f"| {r['target']} | {r['sample_id']} | "
                f"{r['motif_desired_pct']} | {r['alignment_desired_pct']} | "
                f"**{r['difference_pct']}** |"
            )
        report_lines += [
            "",
            "## Interpretation",
            "",
            "Discrepancies between motif-based and alignment-based estimates can arise from:",
            "- Reads with desired motif that also carry indels elsewhere in the read",
            "- Reads where alignment identity falls below the threshold, reducing counted reads",
            "- Differences in read subsampling (--max-reads in 07_alignment_edit_classifier.py)",
            "- Alignment using a short motif reference (not full amplicon) — "
              "indels outside the motif window are not detected",
            "",
            "The motif-count result should be used for all validated biological conclusions.",
        ]
    else:
        report_lines.append("No discrepancies above threshold detected.")

    (out_dir / "discrepancy_report.md").write_text("\n".join(report_lines))

    # Scatter plot (optional)
    if HAS_MPL:
        valid = [(float(r["motif_desired_pct"]),
                  float(r["alignment_desired_pct"]),
                  r["target"], r.get("editor", ""))
                 for r in compare_rows
                 if r["alignment_desired_pct"] not in ("N/A", "")]
        if valid:
            xs, ys, targets, editors = zip(*valid)
            EDITOR_COLORS = {"ePPEmax": "#2166ac", "PE6c": "#d6604d",
                             "WT": "#888888"}
            colors = [EDITOR_COLORS.get(e, "#444") for e in editors]

            fig, ax = plt.subplots(figsize=(5, 5))
            lim = max(max(xs), max(ys)) * 1.1
            ax.plot([0, lim], [0, lim], "k--", linewidth=0.8, label="y = x")
            ax.scatter(xs, ys, c=colors, s=45, edgecolors="black",
                       linewidths=0.5, zorder=3)
            ax.set_xlabel("Motif-count desired % (validated)")
            ax.set_ylabel("Alignment desired % (exploratory, NOT validated)")
            ax.set_title("Motif-count vs alignment-based desired edit estimate\n"
                         "(alignment = NOT validated; motif = gold standard)")
            ax.set_xlim(0, lim); ax.set_ylim(0, lim)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            from matplotlib.patches import Patch
            handles = [Patch(facecolor=c, edgecolor="black", label=e)
                       for e, c in EDITOR_COLORS.items()
                       if e in editors]
            ax.legend(handles=handles, fontsize=8)
            fig.tight_layout()

            for fmt in formats:
                subdir = out_dir / fmt
                subdir.mkdir(exist_ok=True)
                fig.savefig(subdir / f"motif_vs_alignment_scatter.{fmt}",
                            dpi=300, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            print(f"[10_compare] Saved scatter plot in "
                  f"{', '.join(formats)} formats")

    print(f"[10_compare] Outputs in {out_dir}")


if __name__ == "__main__":
    main()
