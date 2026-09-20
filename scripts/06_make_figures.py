#!/usr/bin/env python3
"""
06_make_figures.py  --  Publication-ready figure generation.

*** NEW CONVENIENCE CODE — not part of the validated motif-count analysis. ***

Generates clean, reproducible figures from regenerated pipeline tables.
All biological conclusions come from the validated tables; this script
only visualises them.

Outputs (all under <output>/06_figures/):
  pdf/           vector PDFs
  svg/           vector SVGs
  png/           raster PNGs at --dpi (default 300)
  source_data/   source CSV for every figure
  figure_legends.md  draft figure legends

Seven figures:
  fig1  Desired editing efficiency — main summary (bar + dot)
  fig2  ePPEmax vs PE6c comparison — ALSW (connected dot plot)
  fig3  dpi / time-point effect (grouped bar + dot)
  fig4  WT background / false-positive check
  fig5  Primer-dimer / short-read QC  [informational]
  fig6  Read-pair motif composition (stacked bar)
  fig7  Allele-specific editing — ALSP chr01 vs chr11  [if available]

Style: matplotlib only, no seaborn, white background, open top/right spines.
"""

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils

PIPELINE_ROOT    = Path(__file__).resolve().parent.parent
DEFAULT_METADATA = PIPELINE_ROOT / "config" / "sample_metadata.csv"

# ---------------------------------------------------------------------------
# Optional imports
# ---------------------------------------------------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.lines  as mlines
    import numpy as np
    import pandas as pd
    HAS_DEPS = True
except ImportError as _ie:
    HAS_DEPS = False
    _IMPORT_ERROR = str(_ie)


# ---------------------------------------------------------------------------
# Style and colour constants
# ---------------------------------------------------------------------------

EDITOR_ORDER  = ["WT", "ePPEmax", "PE6c"]
EDITOR_COLORS = {
    "WT":      "#888888",
    "ePPEmax": "#2166ac",
    "PE6c":    "#d6604d",
}
DPI_MARKER = {3: "o", 6: "s"}   # circle = 3 dpi, square = 6 dpi
ALLELE_COLORS = {
    "chr01": "#1a9641",
    "chr11": "#a6d96a",
}
TARGET_DISPLAY = {
    "ALSW":              "ALSW",
    "ALSP_core":         "ALSP (core)",
    "EPSPS_full":        "EPSPS (full edit)",
    "EPSPS_partial1_front": "EPSPS (partial — front)",
    "EPSPS_partial2_back":  "EPSPS (partial — back)",
}

STYLE = {
    "font.family":          "sans-serif",
    "font.size":            10,
    "axes.titlesize":       11,
    "axes.titleweight":     "normal",
    "axes.labelsize":       10,
    "xtick.labelsize":      9,
    "ytick.labelsize":      9,
    "legend.fontsize":      9,
    "legend.frameon":       False,
    "figure.facecolor":     "white",
    "axes.facecolor":       "white",
    "axes.edgecolor":       "black",
    "axes.linewidth":       0.8,
    "axes.grid":            False,
    "axes.spines.top":      False,
    "axes.spines.right":    False,
    "xtick.direction":      "out",
    "ytick.direction":      "out",
    "xtick.major.size":     3,
    "ytick.major.size":     3,
    "xtick.major.width":    0.8,
    "ytick.major.width":    0.8,
    "lines.linewidth":      1.2,
    "patch.linewidth":      0.8,
    "savefig.facecolor":    "white",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_pass_count = 0
_skip_count = 0
_fail_count = 0
_legends: list = []          # accumulated (stem, legend_text) pairs


def _log_pass(stem):
    global _pass_count
    _pass_count += 1
    print(f"[06_figures] SAVED  {stem}")


def _log_skip(stem, reason):
    global _skip_count
    _skip_count += 1
    print(f"[06_figures] SKIP   {stem}: {reason}")


def _log_fail(stem, exc):
    global _fail_count
    _fail_count += 1
    print(f"[06_figures] ERROR  {stem}: {exc}", file=sys.stderr)


def save_fig(fig, stem, out_dir, formats, dpi_val):
    """Save figure in all requested formats; close figure."""
    for fmt in formats:
        subdir = out_dir / fmt
        subdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(subdir / f"{stem}.{fmt}", dpi=dpi_val,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def export_csv(df, stem, out_dir):
    """Write source-data CSV; return path."""
    src_dir = out_dir / "source_data"
    src_dir.mkdir(parents=True, exist_ok=True)
    path = src_dir / f"{stem}.csv"
    df.to_csv(path, index=False)
    return path


def mean_sem(vals):
    """Return (mean, sem) for an array-like; sem=0 if n<2."""
    a = np.asarray([v for v in vals if not math.isnan(float(v))], dtype=float)
    if len(a) == 0:
        return float("nan"), 0.0
    m = float(np.mean(a))
    s = float(np.std(a, ddof=1) / math.sqrt(len(a))) if len(a) > 1 else 0.0
    return m, s


def _editor_x(editors_present, n_groups):
    """Return {editor: x_position} dict with spacing."""
    return {e: i + 1 for i, e in enumerate(editors_present)}


# ---------------------------------------------------------------------------
# Figure 1 — Main editing efficiency
# ---------------------------------------------------------------------------

def fig1_editing_efficiency_main(data, out_dir, formats, dpi_val):
    stem        = "fig1_editing_efficiency_main"
    df_final    = data.get("final_summary")
    if df_final is None or df_final.empty:
        _log_skip(stem, "final_sample_summary.csv not available")
        return

    MAIN_TARGETS = ["ALSW", "ALSP_core", "EPSPS_full"]
    df = df_final[df_final["target"].isin(MAIN_TARGETS)].copy()
    if df.empty:
        _log_skip(stem, "no rows for main targets in final_sample_summary")
        return

    df["pct"]    = pd.to_numeric(df["desired_pct_among_informative"], errors="coerce")
    df["dpi_n"]  = pd.to_numeric(df["dpi"], errors="coerce")

    targets_present = [t for t in MAIN_TARGETS if t in df["target"].values]
    n = len(targets_present)

    fig, axes = plt.subplots(1, n, figsize=(3.8 * n, 4.5), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, tgt in zip(axes, targets_present):
        sub = df[df["target"] == tgt]
        editors = [e for e in EDITOR_ORDER if e in sub["editor"].values]

        for xi, ed in enumerate(editors, start=1):
            esub  = sub[sub["editor"] == ed]
            vals  = esub["pct"].dropna().values
            m, s  = mean_sem(vals)

            # Bar for group mean
            ax.bar(xi, m, width=0.55,
                   color=EDITOR_COLORS.get(ed, "#444"),
                   alpha=0.40, edgecolor="black", linewidth=0.8, zorder=2)

            # SEM error bar
            if s > 0:
                ax.errorbar(xi, m, yerr=s, fmt="none",
                            color="black", capsize=3, linewidth=0.9, zorder=3)

            # Individual sample dots (jitter slightly by dpi)
            jitter_map = {3: -0.10, 6: +0.10}
            for _, row in esub.iterrows():
                dpi_n   = int(row["dpi_n"]) if not math.isnan(row["dpi_n"]) else 0
                jitter  = jitter_map.get(dpi_n, 0.0)
                marker  = DPI_MARKER.get(dpi_n, "D")
                ax.scatter(xi + jitter, row["pct"],
                           color=EDITOR_COLORS.get(ed, "#444"),
                           s=40, zorder=5, marker=marker,
                           edgecolors="black", linewidths=0.5)

        ax.set_xticks(range(1, len(editors) + 1))
        ax.set_xticklabels(editors)
        ax.set_xlim(0.3, len(editors) + 0.7)
        ax.set_ylim(bottom=0)
        ax.set_ylabel("Desired motif frequency (%)")
        ax.set_title(TARGET_DISPLAY.get(tgt, tgt))
        ax.tick_params(axis="x", length=0)

    # Shared legend for dpi markers
    legend_handles = [
        mlines.Line2D([], [], marker="o", color="grey", linestyle="None",
                      markersize=6, label="3 dpi"),
        mlines.Line2D([], [], marker="s", color="grey", linestyle="None",
                      markersize=6, label="6 dpi"),
    ]
    axes[-1].legend(handles=legend_handles, loc="upper right")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.suptitle("Prime editing efficiency by editor and target\n"
                 "(bars = mean of dpi timepoints; dots = individual time points)",
                 fontsize=10, y=1.00)

    try:
        save_fig(fig, stem, out_dir, formats, dpi_val)
        src_cols = ["target", "sample_id", "editor", "dpi",
                    "desired_motif_hits", "informative_reads", "desired_pct_among_informative"]
        export_csv(df[[c for c in src_cols if c in df.columns]].copy(),
                   stem, out_dir)
        _log_pass(stem)
        _legends.append((stem, (
            "Figure 1. Prime editing efficiency by editor and target. "
            "Bar height indicates the mean desired-motif frequency across time points "
            "(3 dpi and 6 dpi). Error bars represent standard error of the mean. "
            "Individual time-point values are overlaid as dots "
            "(circles = 3 dpi; squares = 6 dpi). "
            "Editing frequency is expressed as the number of reads carrying the "
            "desired motif divided by the number of informative reads "
            "(reads carrying either the WT or desired motif). "
            "Targets shown: ALSW (ALS wheat site W), ALSP (ALS wheat site P, core motif), "
            "EPSPS (EPSPS full intended edit). "
            "WT = wild-type control (no prime editor). "
            "Indel and imprecise editing categories are not shown because "
            "they require alignment-based classification not performed in this analysis."
        )))
    except Exception as e:
        _log_fail(stem, e)


# ---------------------------------------------------------------------------
# Figure 2 — ePPEmax vs PE6c comparison (ALSW)
# ---------------------------------------------------------------------------

def fig2_editor_comparison(data, out_dir, formats, dpi_val):
    stem     = "fig2_editor_comparison_ALSW"
    df_final = data.get("final_summary")
    if df_final is None or df_final.empty:
        _log_skip(stem, "final_sample_summary.csv not available")
        return

    df = df_final[
        (df_final["target"] == "ALSW") &
        (df_final["editor"].isin(["ePPEmax", "PE6c"]))
    ].copy()
    if df.empty:
        _log_skip(stem, "no ePPEmax/PE6c ALSW rows")
        return

    # WT reference line
    df_wt = df_final[
        (df_final["target"] == "ALSW") &
        (df_final["editor"] == "WT")
    ].copy()

    df["pct"]   = pd.to_numeric(df["desired_pct_among_informative"], errors="coerce")
    df["dpi_n"] = pd.to_numeric(df["dpi"], errors="coerce")
    dpi_vals    = sorted(df["dpi_n"].dropna().unique())

    fig, ax = plt.subplots(figsize=(5.0, 4.2))

    # WT horizontal reference band
    if not df_wt.empty:
        df_wt["pct"]   = pd.to_numeric(df_wt["desired_pct_among_informative"], errors="coerce")
        wt_mean        = df_wt["pct"].mean()
        ax.axhline(wt_mean, color=EDITOR_COLORS["WT"], linestyle="--",
                   linewidth=1.0, label=f"WT mean ({wt_mean:.2f}%)", zorder=1)

    # One line + dots per editor
    x_pos = {d: i + 1 for i, d in enumerate(dpi_vals)}
    for ed in ["ePPEmax", "PE6c"]:
        esub    = df[df["editor"] == ed].sort_values("dpi_n")
        xs      = [x_pos[d] for d in esub["dpi_n"] if d in x_pos]
        ys      = list(esub["pct"])
        color   = EDITOR_COLORS[ed]
        ax.plot(xs, ys, color=color, linewidth=1.4, zorder=3)
        for xi, yi, dpi_n in zip(xs, ys, esub["dpi_n"]):
            ax.scatter(xi, yi, color=color, s=55, zorder=4,
                       marker=DPI_MARKER.get(int(dpi_n), "D"),
                       edgecolors="black", linewidths=0.5,
                       label=f"{ed} — {int(dpi_n)} dpi")

    ax.set_xticks(list(x_pos.values()))
    ax.set_xticklabels([f"{int(d)} dpi" for d in x_pos])
    ax.set_xlim(0.5, len(dpi_vals) + 0.5)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("Desired motif frequency (%)")
    ax.set_xlabel("Days post-inoculation")
    ax.set_title("ALSW: ePPEmax and PE6c editing efficiency over time")

    # Deduplicated legend
    handles, labels = ax.get_legend_handles_labels()
    seen = {}
    for h, l in zip(handles, labels):
        key = l.split("—")[0].strip()
        if key not in seen:
            seen[key] = (h, l)
    leg_h = [v[0] for v in seen.values()]
    leg_l = [v[1] for v in seen.values()]
    ax.legend(leg_h, leg_l, loc="upper left")

    fig.tight_layout()

    try:
        save_fig(fig, stem, out_dir, formats, dpi_val)
        src_cols = ["target", "sample_id", "editor", "dpi",
                    "desired_motif_hits", "informative_reads", "desired_pct_among_informative"]
        export_csv(df[[c for c in src_cols if c in df.columns]].copy(),
                   stem, out_dir)
        _log_pass(stem)
        _legends.append((stem, (
            "Figure 2. Comparison of ePPEmax and PE6c prime editing efficiency "
            "at the ALSW target over time. "
            "Each data point represents one pooled sample at the indicated time point "
            "(days post-inoculation, dpi). "
            "Lines connect the 3 dpi and 6 dpi values for each editor. "
            "The dashed horizontal line indicates the mean desired-motif frequency "
            "in the wild-type control (WT; no prime editor). "
            "Editing frequency is expressed as desired-motif reads divided by "
            "informative reads (WT + desired motif)."
        )))
    except Exception as e:
        _log_fail(stem, e)


# ---------------------------------------------------------------------------
# Figure 3 — dpi / time-point effect
# ---------------------------------------------------------------------------

def fig3_dpi_effect(data, out_dir, formats, dpi_val):
    stem     = "fig3_dpi_effect_ALSW"
    df_final = data.get("final_summary")
    if df_final is None or df_final.empty:
        _log_skip(stem, "final_sample_summary.csv not available")
        return

    df = df_final[df_final["target"] == "ALSW"].copy()
    if df.empty:
        _log_skip(stem, "no ALSW rows")
        return

    df["pct"]   = pd.to_numeric(df["desired_pct_among_informative"], errors="coerce")
    df["dpi_n"] = pd.to_numeric(df["dpi"], errors="coerce")
    dpi_vals    = sorted(df["dpi_n"].dropna().unique())

    if len(dpi_vals) < 2:
        _log_skip(stem, f"only {len(dpi_vals)} dpi value(s); need ≥2")
        return

    editors     = [e for e in EDITOR_ORDER if e in df["editor"].values]
    n_ed        = len(editors)
    bar_w       = 0.25
    x_base      = np.arange(1, n_ed + 1)
    offsets     = np.linspace(-bar_w * (len(dpi_vals) - 1) / 2,
                               bar_w * (len(dpi_vals) - 1) / 2,
                               len(dpi_vals))

    DPI_HATCHES = {0: "", 1: "///"}

    fig, ax = plt.subplots(figsize=(5.5, 4.2))

    for di, (dv, offset) in enumerate(zip(dpi_vals, offsets)):
        dsub = df[df["dpi_n"] == dv]
        ys   = []
        for ed in editors:
            row = dsub[dsub["editor"] == ed]
            ys.append(float(row["pct"].values[0]) if not row.empty else 0.0)

        bars = ax.bar(x_base + offset, ys,
                      width=bar_w,
                      color=[EDITOR_COLORS.get(e, "#444") for e in editors],
                      edgecolor="black", linewidth=0.7,
                      hatch=DPI_HATCHES.get(di, ""),
                      alpha=0.75,
                      label=f"{int(dv)} dpi",
                      zorder=2)

        # Overlay individual dots
        for xi, yi in zip(x_base + offset, ys):
            ax.scatter(xi, yi, s=30, color="black", zorder=4,
                       edgecolors="black", linewidths=0.5,
                       marker=DPI_MARKER.get(int(dv), "D"))

    ax.set_xticks(x_base)
    ax.set_xticklabels(editors)
    ax.set_xlim(0.4, n_ed + 0.6)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("Desired motif frequency (%)")
    ax.set_xlabel("Editor")
    ax.set_title("ALSW editing efficiency at 3 dpi and 6 dpi\n"
                 "(hatched = 3 dpi, solid = 6 dpi)")
    ax.tick_params(axis="x", length=0)

    # Legend: dpi times
    dpi_legend = [
        mpatches.Patch(facecolor="grey", hatch=DPI_HATCHES[i],
                       edgecolor="black", label=f"{int(dv)} dpi")
        for i, dv in enumerate(dpi_vals)
    ]
    ax.legend(handles=dpi_legend, loc="upper right")

    fig.tight_layout()

    try:
        save_fig(fig, stem, out_dir, formats, dpi_val)
        src_cols = ["target", "sample_id", "editor", "dpi",
                    "desired_motif_hits", "informative_reads", "desired_pct_among_informative"]
        export_csv(df[[c for c in src_cols if c in df.columns]].copy(),
                   stem, out_dir)
        _log_pass(stem)
        _legends.append((stem, (
            "Figure 3. Effect of time point (days post-inoculation, dpi) on ALSW "
            "prime editing efficiency. "
            "Hatched bars represent 3 dpi; solid bars represent 6 dpi. "
            "Each bar represents one pooled sample. "
            "Editing frequency is expressed as desired-motif reads divided by "
            "informative reads (WT + desired motif). "
            "WT = wild-type control (no prime editor)."
        )))
    except Exception as e:
        _log_fail(stem, e)


# ---------------------------------------------------------------------------
# Figure 4 — WT background / false-positive check
# ---------------------------------------------------------------------------

def fig4_wt_background(data, out_dir, formats, dpi_val):
    stem   = "fig4_wt_background"
    df_wt  = data.get("wt_background")
    if df_wt is None or df_wt.empty:
        _log_skip(stem, "wt_background_summary.csv not available")
        return

    df     = df_wt.copy()
    df["pct"] = pd.to_numeric(df["background_pct"], errors="coerce")
    df["dpi"] = pd.to_numeric(df["dpi"], errors="coerce")

    # Only show main targets to avoid over-crowding
    SHOW  = ["ALSW", "ALSP_core", "EPSPS_full"]
    df    = df[df["target"].isin(SHOW)]
    if df.empty:
        _log_skip(stem, "no rows for main targets")
        return

    df = df.sort_values(["target", "dpi"])
    labels = [f"{r['target']}\n({int(r['dpi'])} dpi)" for _, r in df.iterrows()]
    xs     = range(len(df))
    colors = [EDITOR_COLORS["WT"]] * len(df)

    fig, ax = plt.subplots(figsize=(max(6.0, len(df) * 0.65), 4.2))

    ax.bar(xs, df["pct"], color=colors,
           edgecolor="black", linewidth=0.8, width=0.6)
    ax.scatter(xs, df["pct"], s=35, color="black", zorder=4)

    ax.set_xticks(list(xs))
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_ylabel("WT background desired-motif frequency (%)")
    ax.set_title("WT control: desired-motif background signal by target\n"
                 "(motif-level background; not alignment-confirmed editing)")
    ax.set_ylim(bottom=0)

    # Annotate the ALSW bars since they are non-zero
    for xi, (_, row) in zip(xs, df.iterrows()):
        if row["pct"] > 0:
            ax.text(xi, row["pct"] + 0.02 * ax.get_ylim()[1],
                    f"{row['pct']:.2f}%", ha="center", va="bottom",
                    fontsize=8)

    fig.tight_layout()

    try:
        save_fig(fig, stem, out_dir, formats, dpi_val)
        src_cols = ["target", "sample_id", "dpi",
                    "wt_motif_hits", "desired_motif_hits",
                    "informative_reads", "background_pct"]
        export_csv(df[[c for c in src_cols if c in df.columns]].copy(),
                   stem, out_dir)
        _log_pass(stem)
        _legends.append((stem, (
            "Figure 4. Desired-motif background signal in wild-type (WT) control "
            "samples (no prime editor). "
            "Background is expressed as the number of reads matching the "
            "desired motif divided by informative reads (WT + desired motif) "
            "in the WT sample. "
            "Non-zero ALSW background (~0.37–0.57%) may reflect PCR-introduced "
            "substitutions, sequencing errors, or biological variation at this "
            "site; it should be subtracted from editor samples for conservative "
            "efficiency estimates. "
            "ALSP and EPSPS show no detectable background at this motif resolution. "
            "This is a motif-level signal; alignment-based variant calling was "
            "not performed."
        )))
    except Exception as e:
        _log_fail(stem, e)


# ---------------------------------------------------------------------------
# Figure 5 — Primer-dimer / short-read QC
# ---------------------------------------------------------------------------

def fig5_primer_dimer_qc(data, out_dir, formats, dpi_val):
    stem    = "fig5_primer_dimer_qc"
    df_dim  = data.get("dimer_qc")
    if df_dim is None or df_dim.empty:
        _log_skip(stem, "primer_dimer_qc_summary.csv not available")
        return

    df = df_dim.copy()

    # Rename sample column for consistency (02_qc uses 'sample', not 'sample_id')
    if "sample" in df.columns and "sample_id" not in df.columns:
        df = df.rename(columns={"sample": "sample_id"})

    df["short_pct"]     = pd.to_numeric(df.get("short_fraction_pct",   df.get("short_pct", 0)), errors="coerce")
    df["amplicon_pct"]  = pd.to_numeric(df.get("amplicon_fraction_pct", pd.Series(dtype=float)), errors="coerce")

    # Exclude undetermined
    df = df[df["sample_id"] != "undetermined"].copy()
    df = df.sort_values("short_pct", ascending=False)

    # Add editor metadata if available
    meta = data.get("metadata")
    if meta is not None:
        df = df.merge(meta[["sample_id", "editor", "target"]],
                      on="sample_id", how="left")

    colors = [EDITOR_COLORS.get(e, "#aaaaaa") for e in df.get("editor", [""] * len(df))]

    fig, ax = plt.subplots(figsize=(max(8.0, len(df) * 0.55), 4.2))

    xs = range(len(df))
    ax.bar(xs, df["short_pct"], color=colors,
           edgecolor="black", linewidth=0.7, width=0.7, zorder=2)
    ax.axhline(20, color="red", linestyle="--", linewidth=1.0,
               label="20% reference line", zorder=3)

    # Useful-read fraction as a secondary reference
    if df["amplicon_pct"].notna().any():
        ax.plot(list(xs), list(df["amplicon_pct"]), "k^",
                markersize=5, zorder=4, label="Amplicon fraction (%)")

    ax.set_xticks(list(xs))
    ax.set_xticklabels(df["sample_id"], rotation=55, ha="right", fontsize=8)
    ax.set_ylabel("Short-read / primer-dimer fraction (%)")
    ax.set_ylim(0, 105)
    ax.set_title("Short-read / primer-dimer fraction by sample  [QC — reads not filtered]")

    # Editor colour legend
    ed_handles = [mpatches.Patch(facecolor=EDITOR_COLORS[e], edgecolor="black",
                                 label=e)
                  for e in EDITOR_ORDER if e in df.get("editor", pd.Series()).values]
    ed_handles.append(mlines.Line2D([], [], color="red", linestyle="--",
                                    label="20% reference"))
    ax.legend(handles=ed_handles, loc="upper right", fontsize=8)

    fig.tight_layout()

    try:
        save_fig(fig, stem, out_dir, formats, dpi_val)
        src_cols = ["sample_id"]
        for c in ["editor", "target", "total_read_pairs",
                  "short_reads_lt_threshold", "short_read_threshold_bp",
                  "short_fraction_pct", "amplicon_fraction_pct"]:
            if c in df.columns:
                src_cols.append(c)
        export_csv(df[src_cols].copy(), stem, out_dir)
        _log_pass(stem)
        _legends.append((stem, (
            "Figure 5. Short-read / primer-dimer fraction by sample. "
            "Bars show the fraction of reads shorter than the minimum amplicon "
            "length threshold (100 bp post-trim) in each demultiplexed sample. "
            "The red dashed line marks 20%. "
            "Triangles indicate the amplicon-length fraction (1 − short fraction). "
            "This figure is informational QC only: reads were NOT filtered before "
            "motif counting. The validated editing percentages in Figures 1–3 are "
            "based on informative reads (WT + desired motif), not total read pairs, "
            "so the high short-read fraction does not invalidate those results. "
            "Colours indicate editor; bars from all three editor groups are shown."
        )))
    except Exception as e:
        _log_fail(stem, e)


# ---------------------------------------------------------------------------
# Figure 6 — Read-pair motif composition (stacked bar)
# ---------------------------------------------------------------------------

def fig6_motif_composition(data, out_dir, formats, dpi_val):
    stem     = "fig6_motif_composition"
    df_motif = data.get("motif_summary")
    if df_motif is None or df_motif.empty:
        _log_skip(stem, "motif_count_summary.tsv not available")
        return

    SHOW_TARGETS = ["ALSW", "ALSP_core", "EPSPS_full"]
    df = df_motif[df_motif["target"].isin(SHOW_TARGETS)].copy()
    if df.empty:
        _log_skip(stem, "no rows for selected targets")
        return

    for col in ["total_read_pairs", "wt_motif_hits",
                "desired_motif_hits", "both_hits", "neither_hits"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

    # Fractions relative to total
    tot = df["total_read_pairs"].replace(0, np.nan)
    df["frac_wt"]      = df["wt_motif_hits"]      / tot * 100
    df["frac_desired"] = df["desired_motif_hits"]  / tot * 100
    df["frac_both"]    = df["both_hits"]           / tot * 100
    df["frac_neither"] = df["neither_hits"]        / tot * 100

    n_targets = len(SHOW_TARGETS)
    fig, axes = plt.subplots(1, n_targets, figsize=(4.5 * n_targets, 4.5),
                              sharey=True)
    if n_targets == 1:
        axes = [axes]

    STACK_COLORS   = ["#4393c3", "#2166ac", "#fdae61", "#d1d1d1"]
    STACK_LABELS   = ["WT motif", "Desired motif", "Both motifs", "Neither"]
    STACK_FRACS    = ["frac_wt", "frac_desired", "frac_both", "frac_neither"]

    for ax, tgt in zip(axes, SHOW_TARGETS):
        sub = df[df["target"] == tgt].reset_index(drop=True)
        if sub.empty:
            ax.set_visible(False)
            continue

        xs      = range(len(sub))
        bottoms = [0.0] * len(sub)

        for frac_col, color, label in zip(STACK_FRACS, STACK_COLORS, STACK_LABELS):
            vals = sub[frac_col].fillna(0).values
            ax.bar(xs, vals, bottom=bottoms,
                   color=color, edgecolor="none", width=0.7,
                   label=label)
            bottoms = [b + v for b, v in zip(bottoms, vals)]

        ax.set_xticks(list(xs))
        ax.set_xticklabels(sub["sample"], rotation=55, ha="right", fontsize=7)
        ax.set_title(TARGET_DISPLAY.get(tgt, tgt))
        ax.set_ylim(0, 105)

    axes[0].set_ylabel("Fraction of total read pairs (%)")

    # Single shared legend on right panel
    handles = [mpatches.Patch(facecolor=c, label=l)
               for c, l in zip(STACK_COLORS, STACK_LABELS)]
    axes[-1].legend(handles=handles, loc="upper right", fontsize=8)

    fig.suptitle("Read-pair motif composition by sample\n"
                 "(indel/imprecise categories not shown — alignment not performed)",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    try:
        save_fig(fig, stem, out_dir, formats, dpi_val)
        src_cols = ["target", "sample", "total_read_pairs",
                    "wt_motif_hits", "desired_motif_hits",
                    "both_hits", "neither_hits",
                    "frac_wt", "frac_desired", "frac_both", "frac_neither"]
        export_csv(df[[c for c in src_cols if c in df.columns]].copy(),
                   stem, out_dir)
        _log_pass(stem)
        _legends.append((stem, (
            "Figure 6. Motif composition of sequenced read pairs by sample and target. "
            "Stacked bars show the fraction of total read pairs carrying the WT motif "
            "(light blue), desired edit motif (dark blue), both motifs simultaneously "
            "(orange; expected to be near zero), or neither motif (grey). "
            "Read pairs without either motif typically represent short/non-specific "
            "amplification products and do not contribute to editing efficiency estimates. "
            "Indel and imprecise prime-editing categories are not shown because "
            "alignment-based variant classification was not performed in this analysis. "
            "Targets shown: ALSW, ALSP core motif, EPSPS full intended edit."
        )))
    except Exception as e:
        _log_fail(stem, e)


# ---------------------------------------------------------------------------
# Figure 7 — Allele-specific editing (ALSP chr01 vs chr11)
# ---------------------------------------------------------------------------

def fig7_allele_specific(data, out_dir, formats, dpi_val):
    stem       = "fig7_allele_specific_ALSP"
    df_allele  = data.get("allele_summary")
    if df_allele is None or df_allele.empty:
        _log_skip(stem, "allele_specific_summary.csv not available")
        return

    df = df_allele[df_allele["allele"].isin(["chr01", "chr11"])].copy()
    if df.empty:
        _log_skip(stem, "no chr01/chr11 allele rows (ALSW not allele-split — "
                        "validated design)")
        return

    df["pct"]   = pd.to_numeric(df["desired_pct"], errors="coerce")
    df["dpi_n"] = pd.to_numeric(df["dpi"], errors="coerce")

    alleles  = ["chr01", "chr11"]
    editors  = [e for e in EDITOR_ORDER if e in df["editor"].values]
    bar_w    = 0.30
    x_base   = np.arange(1, len(editors) + 1)
    offsets  = [-bar_w / 2, +bar_w / 2]

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.2), sharey=True)

    dpi_vals = sorted(df["dpi_n"].dropna().unique())
    for ax, dv in zip(axes, dpi_vals[:2]):
        sub = df[df["dpi_n"] == dv]
        for al, offset in zip(alleles, offsets):
            asub = sub[sub["allele"] == al]
            ys   = []
            for ed in editors:
                row = asub[asub["editor"] == ed]
                ys.append(float(row["pct"].values[0]) if not row.empty else 0.0)

            ax.bar(x_base + offset, ys,
                   width=bar_w,
                   color=ALLELE_COLORS[al],
                   edgecolor="black", linewidth=0.7,
                   label=al, alpha=0.85)

        ax.set_xticks(x_base)
        ax.set_xticklabels(editors)
        ax.set_xlim(0.4, len(editors) + 0.6)
        ax.set_ylim(bottom=0)
        ax.set_title(f"{int(dv)} dpi")
        ax.tick_params(axis="x", length=0)

    axes[0].set_ylabel("Desired motif frequency (%)")
    fig.suptitle("ALSP allele-specific prime editing efficiency "
                 "(chr01 vs chr11)\n"
                 "(allele discrimination based on flanking-SNP motifs)",
                 fontsize=10)

    allele_handles = [mpatches.Patch(facecolor=ALLELE_COLORS[a],
                                     edgecolor="black", label=a)
                      for a in alleles]
    axes[-1].legend(handles=allele_handles, loc="upper right")

    fig.tight_layout(rect=[0, 0, 1, 0.90])

    try:
        save_fig(fig, stem, out_dir, formats, dpi_val)
        src_cols = ["target", "allele", "sample_id", "editor", "dpi",
                    "wt_hits", "desired_hits", "informative", "desired_pct"]
        export_csv(df[[c for c in src_cols if c in df.columns]].copy(),
                   stem, out_dir)
        _log_pass(stem)
        _legends.append((stem, (
            "Figure 7. Allele-specific prime editing efficiency at the ALSP target "
            "(chromosome 1 vs chromosome 11 homeologs). "
            "Allele discrimination is based on flanking single-nucleotide "
            "polymorphisms (SNPs) in the motif sequence: "
            "TTGTTGCTATAACT (chr01-specific) vs TTGTTGCTATAACCC (chr11-specific). "
            "Bars show the fraction of allele-assigned reads carrying the desired "
            "motif for each editor and time point. "
            "ALSW is not shown because no validated allele-discriminating motifs "
            "exist for that target in this analysis. "
            "ALSP editing values are generally low; interpret with caution "
            "(low desired-read counts)."
        )))
    except Exception as e:
        _log_fail(stem, e)


# ---------------------------------------------------------------------------
# Figure 8 — Allele-aware precise edit (optional; needs 11_allele_aware outputs)
# ---------------------------------------------------------------------------

def fig8_allele_aware_precise_edit(data, out_dir, formats, dpi_val):
    """
    Fig 8: Allele-aware precise edit quantification (chr01 vs chr11).
    Only generated if 11_allele_aware/allele_aware_precise_edit_summary.csv exists.
    Reads allele-aware data from the allele_aware key in data dict.

    *** NOT VALIDATED — exploratory only ***
    """
    stem       = "fig8_allele_aware_precise_edit"
    df_aa      = data.get("allele_aware_summary")

    if df_aa is None or df_aa.empty:
        _log_skip(stem, "allele_aware_precise_edit_summary.csv not found — "
                        "run with --run-allele-aware to generate")
        return

    # Filter to actual allele rows (not ambiguous/low_information)
    df = df_aa[df_aa["allele"].isin(["chr01", "chr11"])].copy()
    if df.empty:
        _log_skip(stem, "no chr01/chr11 rows in allele-aware summary")
        return

    df["pct"]   = pd.to_numeric(df["precise_desired_pct_of_allele_assigned"],
                                 errors="coerce")
    df["dpi_n"] = pd.to_numeric(df["dpi"], errors="coerce")

    alleles  = ["chr01", "chr11"]
    editors  = [e for e in ["WT", "ePPEmax", "PE6c"] if e in df["editor"].values]
    dpi_vals = sorted(df["dpi_n"].dropna().unique())
    n_dpi    = min(len(dpi_vals), 2)

    if n_dpi == 0:
        _log_skip(stem, "no dpi values in allele-aware data")
        return

    bar_w   = 0.30
    x_base  = np.arange(1, len(editors) + 1)
    offsets = [-bar_w / 2, +bar_w / 2]

    fig, axes = plt.subplots(1, n_dpi, figsize=(5 * n_dpi, 4.2), sharey=True)
    if n_dpi == 1:
        axes = [axes]

    for ax, dv in zip(axes, dpi_vals[:2]):
        sub = df[df["dpi_n"] == dv]
        for al, offset in zip(alleles, offsets):
            asub = sub[sub["allele"] == al]
            ys = []
            for ed in editors:
                row = asub[asub["editor"] == ed]
                ys.append(float(row["pct"].values[0]) if not row.empty else 0.0)
            ax.bar(x_base + offset, ys,
                   width=bar_w,
                   color=ALLELE_COLORS.get(al, "#888"),
                   edgecolor="black", linewidth=0.7,
                   label=al, alpha=0.85)

        ax.set_xticks(x_base)
        ax.set_xticklabels(editors)
        ax.set_xlim(0.4, len(editors) + 0.6)
        ax.set_ylim(bottom=0)
        ax.set_title(f"{int(dv)} dpi")
        ax.tick_params(axis="x", length=0)

    axes[0].set_ylabel("Precise desired edit % of allele-assigned reads\n"
                       "(allele-aware, NOT validated)")
    fig.suptitle("ALSP allele-aware precise edit quantification"
                 " (chr01 vs chr11)\n"
                 "*** EXPLORATORY — not validated against controls ***",
                 fontsize=10)

    allele_handles = [
        mpatches.Patch(facecolor=ALLELE_COLORS.get(a, "#888"),
                       edgecolor="black", label=a)
        for a in alleles
    ]
    axes[-1].legend(handles=allele_handles, loc="upper right")
    fig.tight_layout(rect=[0, 0, 1, 0.88])

    try:
        save_fig(fig, stem, out_dir, formats, dpi_val)
        src_cols = [c for c in ["target", "allele", "sample_id", "editor", "dpi",
                                 "allele_assigned_reads", "precise_desired_reads",
                                 "precise_desired_pct_of_allele_assigned"]
                    if c in df.columns]
        export_csv(df[src_cols].copy(), stem, out_dir)
        _log_pass(stem)
        _legends.append((stem, (
            "Figure 8. Allele-aware precise edit quantification at the ALSP target. "
            "Allele assignment is based on a single allele-discriminating SNP at position 13 "
            "of the 29 bp context motif (chr01: T; chr11: C). "
            "The edit position (position 23; C→T) is excluded from allele assignment. "
            "Bars show the percentage of allele-assigned reads carrying the precise "
            "desired edit, split by allele (chr01 = dark green; chr11 = light green). "
            "*** IMPORTANT: These data are exploratory and have NOT been validated against "
            "known controls. Do not use in manuscripts without independent validation. "
            "The validated editing percentages are in Figure 7 (motif-count method). ***"
        )))
    except Exception as e:
        _log_fail(stem, e)


# ---------------------------------------------------------------------------
# Figure legends markdown
# ---------------------------------------------------------------------------

def write_figure_legends(out_dir):
    path = out_dir / "figure_legends.md"
    lines = [
        "# Figure Legends",
        "",
        "> These legends were generated automatically by `06_make_figures.py`.",
        "> Review and edit before submission.",
        "> Caveats about the analysis method are included where relevant.",
        "",
    ]
    for stem, legend in _legends:
        lines.append(f"## {stem}")
        lines.append("")
        lines.append(legend)
        lines.append("")
    path.write_text("\n".join(lines))
    print(f"[06_figures] SAVED  figure_legends.md")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", required=True,
                   help="Root output directory containing 02_qc/, 03_allele/, "
                        "04_edit/, 05_summary/")
    p.add_argument("--output", required=True,
                   help="Root output directory (figures go into <output>/06_figures/)")
    p.add_argument("--metadata", default=str(DEFAULT_METADATA),
                   help="sample_metadata.csv path")
    p.add_argument("--input-summary-dir", default=None, dest="summary_dir",
                   help="Override path to summary tables directory "
                        "[default: <input>/05_summary/]")
    p.add_argument("--format", default="pdf,png,svg",
                   help="Comma-separated output formats [default: pdf,png,svg]")
    p.add_argument("--dpi", type=int, default=300,
                   help="Raster DPI for PNG (and any raster format) [default: 300]")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    if not HAS_DEPS:
        print(f"[06_figures] ERROR: missing dependencies: {_IMPORT_ERROR}",
              file=sys.stderr)
        print("[06_figures]   Run: pip install matplotlib pandas numpy",
              file=sys.stderr)
        sys.exit(1)

    plt.rcParams.update(STYLE)

    formats = [f.strip().lower() for f in args.format.split(",")]
    dpi_val = args.dpi
    in_root = Path(args.input)
    out_dir = Path(args.output) / "06_figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_dir = Path(args.summary_dir) if args.summary_dir else in_root / "05_summary"

    print(f"[06_figures] Input root   : {in_root}")
    print(f"[06_figures] Summary dir  : {summary_dir}")
    print(f"[06_figures] Output dir   : {out_dir}")
    print(f"[06_figures] Formats      : {formats}")
    print(f"[06_figures] PNG DPI      : {dpi_val}")

    # --- Load all data tables ---
    def _csv(path):
        try:
            return pd.read_csv(path)
        except FileNotFoundError:
            return None

    def _tsv(path):
        try:
            return pd.read_csv(path, sep="\t")
        except FileNotFoundError:
            return None

    meta_path = Path(args.metadata)
    data = {
        "metadata":             _csv(meta_path)  if meta_path.exists() else None,
        "final_summary":        _csv(summary_dir / "final_sample_summary.csv"),
        "wt_background":        _csv(summary_dir / "wt_background_summary.csv"),
        "allele_summary":       _csv(summary_dir / "allele_specific_summary.csv"),
        "dimer_qc":             _csv(in_root / "02_qc" / "primer_dimer_qc_summary.csv"),
        "motif_summary":        _tsv(in_root / "04_edit" / "motif_count_summary.tsv"),
        # Optional: allele-aware summary from 11_allele_aware_precise_edit.py
        "allele_aware_summary": _csv(
            in_root / "11_allele_aware" / "allele_aware_precise_edit_summary.csv"
        ),
    }

    missing = [k for k, v in data.items() if v is None]
    if missing:
        print(f"[06_figures] Tables not found (will skip affected figures): "
              f"{missing}")

    # --- Generate figures ---
    fig1_editing_efficiency_main(data, out_dir, formats, dpi_val)
    fig2_editor_comparison(      data, out_dir, formats, dpi_val)
    fig3_dpi_effect(             data, out_dir, formats, dpi_val)
    fig4_wt_background(          data, out_dir, formats, dpi_val)
    fig5_primer_dimer_qc(        data, out_dir, formats, dpi_val)
    fig6_motif_composition(      data, out_dir, formats, dpi_val)
    fig7_allele_specific(        data, out_dir, formats, dpi_val)
    fig8_allele_aware_precise_edit(data, out_dir, formats, dpi_val)

    # --- Legends ---
    write_figure_legends(out_dir)

    # --- Summary ---
    all_files = sorted(out_dir.rglob("*.*"))
    print(f"\n[06_figures] Saved {len(all_files)} files in {out_dir}/")
    for f in all_files:
        print(f"  {f.relative_to(out_dir)}")

    print(f"\n[06_figures] Figures saved:  {_pass_count}")
    print(f"[06_figures] Figures skipped: {_skip_count}")
    if _fail_count:
        print(f"[06_figures] Figures FAILED:  {_fail_count}", file=sys.stderr)


if __name__ == "__main__":
    main()
