#!/usr/bin/env python3
"""
11_allele_aware_precise_edit.py  --  Allele-aware precise edit quantification.

*** OPTIONAL MODULE — NOT PART OF THE VALIDATED MOTIF-COUNT PIPELINE ***
*** OUTPUTS ARE NEW AND REQUIRE INDEPENDENT VALIDATION                 ***

Purpose
-------
For targets with multiple homeologous alleles (e.g. Nicotiana benthamiana
ALS chr01 and chr11), this module assigns each read to a specific allele
using allele-discriminating SNP positions, then independently calls edit
status within each allele.

The validated motif_count_summary.tsv is NEVER modified.

Algorithm
---------
For each read pair:
  1. Build combined search string:
       R1 + N + R2 + N + RC(R1) + N + RC(R2)
  2. Search for allele_common_anchor in the combined string.
  3. The base immediately after the anchor is the allele-discriminating SNP.
     Check this base against the expected bases for chr01/chr11.
  4. The edit base is at anchor_end + SNP_offset + edit_offset_from_snp.
     (Edit position is EXCLUDED from allele assignment.)
  5. Assign allele (chr01 / chr11 / ambiguous / low_information).
  6. Only for allele-assigned reads, call edit status
     (wt / precise_desired / non_desired / ambiguous_edit).

Allele assignment categories
-----------------------------
  chr01             SNP base matches chr01 expected base(s) unambiguously
  chr11             SNP base matches chr11 expected base(s) unambiguously
  ambiguous         SNP evidence tied between alleles (e.g. both match or
                    neither matches uniquely given score margin)
  low_information   Anchor not found in read (SNP site not covered)

Edit status categories (allele-assigned reads only)
-----------------------------------------------------
  wt                Edit base = WT base
  precise_desired   Edit base = desired base
  non_desired       Edit base present but not WT and not desired
  ambiguous_edit    Edit site not covered (read ends before edit position)

Outputs (under <output>/11_allele_aware/):
  allele_aware_precise_edit_summary.csv      ← main result table
  read_level_allele_edit_calls.tsv           ← per-read detail
  allele_assignment_qc.csv                   ← assignment rate QC
  allele_aware_vs_motif_comparison.csv       ← cross-check with motif-count
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

MAX_READS_DEFAULT = 0   # 0 = process all reads

# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------

def _find_anchor_and_base(combined: str, anchor: str, offset: int):
    """
    Find anchor in combined string.
    Returns (base_at_anchor_end_plus_offset, anchor_found_position)
    or (None, -1) if anchor not found.
    """
    pos = combined.find(anchor.upper())
    if pos < 0:
        return None, -1
    target_pos = pos + len(anchor) + offset
    if target_pos >= len(combined) or combined[target_pos] == 'N':
        return None, pos
    return combined[target_pos], pos


def analyze_read(s1: str, s2: str,
                 anchor: str,
                 allele_snps: dict,       # {allele_name: [{"snp_offset": int, "expected_base": str}]}
                 edit_offset_from_snp: int,
                 wt_base: str,
                 desired_base: str,
                 min_snp_score_margin: int = 1,
                 min_informative_snps: int = 1,
                 edit_window_start: int = None,
                 edit_window_end: int = None,
                 ) -> dict:
    """
    Assign allele and call edit for a single read pair.

    The edit base is at position: anchor_end + snp_offset + edit_offset_from_snp
    where snp_offset is always 0 in the current ALSP config.

    Returns a flat dict of classification fields.
    """
    combined = (s1.upper() + "N"
                + s2.upper() + "N"
                + utils.reverse_complement(s1) + "N"
                + utils.reverse_complement(s2))

    # ── Step 1: find anchor ──────────────────────────────────────────────
    anchor_pos = combined.find(anchor.upper())
    if anchor_pos < 0:
        return _low_info_result("anchor_not_found_in_read")

    anchor_end = anchor_pos + len(anchor)   # position of SNP in combined

    # ── Step 2: score each allele against its SNPs ───────────────────────
    allele_scores   = {}
    allele_covered  = {}   # n SNP positions that were readable

    for allele_name, snp_list in allele_snps.items():
        score    = 0
        covered  = 0
        for snp_def in snp_list:
            snp_idx = anchor_end + snp_def["snp_offset"]
            if snp_idx >= len(combined) or combined[snp_idx] == 'N':
                continue   # position not covered
            covered += 1
            obs = combined[snp_idx].upper()
            if obs == snp_def["expected_base"].upper():
                score += 1
        allele_scores[allele_name]  = score
        allele_covered[allele_name] = covered

    max_informative = max(allele_covered.values(), default=0)
    if max_informative < min_informative_snps:
        return _low_info_result("insufficient_snp_coverage")

    # ── Step 3: assign allele from scores ────────────────────────────────
    sorted_alleles = sorted(allele_scores, key=allele_scores.get, reverse=True)
    if not sorted_alleles:
        return _low_info_result("no_alleles_configured")
    best = sorted_alleles[0]
    # A config with a single allele has no runner-up to compare against; the
    # margin is then the winner's own score (indexing [1] raised IndexError).
    second = sorted_alleles[1] if len(sorted_alleles) > 1 else None
    margin = (allele_scores[best] - allele_scores[second]
              if second is not None else allele_scores[best])

    if allele_scores[best] == 0:
        assigned_allele = "ambiguous"
        ambiguous_reason = "no_snp_matches_any_allele"
    elif margin < min_snp_score_margin:
        assigned_allele  = "ambiguous"
        ambiguous_reason = f"margin_{margin}_below_threshold_{min_snp_score_margin}"
    else:
        assigned_allele  = best
        ambiguous_reason = ""

    # SNP base for reporting (first SNP, offset 0)
    first_snp_idx  = anchor_end   # offset 0
    snp_base_obs   = (combined[first_snp_idx].upper()
                      if first_snp_idx < len(combined) and combined[first_snp_idx] != 'N'
                      else "N")
    snp_score_chr01 = allele_scores.get("chr01", 0)
    snp_score_chr11 = allele_scores.get("chr11", 0)
    n_informative   = max_informative

    # ── Step 4: call edit (only meaningful for assigned reads) ───────────
    edit_idx = anchor_end + edit_offset_from_snp
    if edit_idx >= len(combined) or combined[edit_idx] == 'N':
        edit_base_obs = "N"
        edit_status   = "ambiguous_edit"
    else:
        edit_base_obs = combined[edit_idx].upper()
        if edit_base_obs == desired_base.upper():
            edit_status = "precise_desired"
        elif edit_base_obs == wt_base.upper():
            edit_status = "wt"
        else:
            edit_status = "non_desired"

    has_precise = (edit_status == "precise_desired" and
                   assigned_allele not in ("ambiguous", "low_information"))

    return {
        "assigned_allele":           assigned_allele,
        "allele_score_chr01":        snp_score_chr01,
        "allele_score_chr11":        snp_score_chr11,
        "informative_snp_count":     n_informative,
        "snp_base_observed":         snp_base_obs,
        "ambiguous_reason":          ambiguous_reason,
        "edit_position_base_obs":    edit_base_obs,
        "edit_status":               edit_status,
        "has_precise_desired_edit":  has_precise,
    }


def _low_info_result(reason: str) -> dict:
    return {
        "assigned_allele":           "low_information",
        "allele_score_chr01":        0,
        "allele_score_chr11":        0,
        "informative_snp_count":     0,
        "snp_base_observed":         "N",
        "ambiguous_reason":          reason,
        "edit_position_base_obs":    "N",
        "edit_status":               "ambiguous_edit",
        "has_precise_desired_edit":  False,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input",    required=True,
                   help="Root output dir (expects <input>/01_demux/demux_clean/)")
    p.add_argument("--output",   required=True,
                   help="Root output dir (outputs to <output>/11_allele_aware/)")
    p.add_argument("--metadata", default=str(DEFAULT_METADATA))
    p.add_argument("--config",   default=str(DEFAULT_CONFIG))
    p.add_argument("--target",   default=None,
                   help="Restrict to one allele-aware target name (e.g. ALSP). "
                        "Default: all configured targets.")
    p.add_argument("--min-informative-snps", type=int, default=1,
                   dest="min_informative_snps",
                   help="Min SNP positions covered to attempt assignment [default: 1]")
    p.add_argument("--min-allele-margin", type=int, default=1,
                   dest="min_allele_margin",
                   help="Min score difference to make an unambiguous allele call "
                        "[default: 1]")
    p.add_argument("--max-reads", type=int, default=MAX_READS_DEFAULT,
                   dest="max_reads",
                   help="Max reads per sample (0 = all) [default: 0 = all]")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args      = parse_args()
    clean_dir = Path(args.input) / "01_demux" / "demux_clean"
    out_dir   = Path(args.output) / "11_allele_aware"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not clean_dir.exists():
        print(f"[11_allele] ERROR: demux_clean not found: {clean_dir}",
              file=sys.stderr)
        sys.exit(1)

    # Load configs
    meta_rows       = utils.read_csv_file(args.metadata)
    sample_to_meta  = {r["sample_id"]: r for r in meta_rows}
    full_cfg        = utils.load_json(args.config)
    aa_targets      = full_cfg.get("allele_aware_targets", [])

    if not aa_targets:
        print("[11_allele] No allele_aware_targets in targets.json — nothing to do.")
        print("[11_allele] Add allele_aware_targets entries to targets.json "
              "to enable this module.")
        return

    # Filter to requested target if specified
    if args.target:
        aa_targets = [t for t in aa_targets if t["target_name"] == args.target]
        if not aa_targets:
            print(f"[11_allele] ERROR: allele-aware target '{args.target}' not found.",
                  file=sys.stderr)
            sys.exit(1)

    print(f"[11_allele] *** NOT VALIDATED — exploratory analysis only ***")

    # Accumulate across all targets for combined output files
    all_read_rows    = []
    all_sample_rows  = []
    all_qc_rows      = []

    for aa_cfg in aa_targets:
        t_name = aa_cfg["target_name"]
        status = aa_cfg.get("status", "unconfigured")

        if status != "configured":
            print(f"[11_allele] SKIP {t_name}: status='{status}' — "
                  f"{aa_cfg.get('notes', '')[:80]}")
            continue

        anchor          = aa_cfg["allele_common_anchor"].upper()
        allele_cfg      = aa_cfg["alleles"]
        edit_def        = aa_cfg["edit"]
        edit_offset     = edit_def["offset_from_snp"]
        wt_base         = edit_def["wt_base"].upper()
        desired_base    = edit_def["desired_base"].upper()
        applicable      = [a.upper() for a in aa_cfg["applicable_to"]]
        edit_win_start  = aa_cfg.get("edit_window", {}).get("start")
        edit_win_end    = aa_cfg.get("edit_window", {}).get("end")

        # Reshape allele SNPs for the algorithm
        allele_snps = {
            a_name: [
                {"snp_offset": s["snp_offset_from_anchor_end"],
                 "expected_base": s["expected_base"]}
                for s in a_data["allele_snps"]
            ]
            for a_name, a_data in allele_cfg.items()
        }
        allele_names = list(allele_snps.keys())

        print(f"[11_allele] Target: {t_name}  "
              f"anchor={anchor}  alleles={allele_names}  "
              f"edit_offset_from_snp={edit_offset}")

        # Iterate samples in metadata order
        for row in meta_rows:
            sid     = row["sample_id"]
            tgene   = row["target"].upper()
            if tgene not in applicable:
                continue

            r1 = clean_dir / f"{sid}_R1.fastq.gz"
            r2 = clean_dir / f"{sid}_R2.fastq.gz"
            if not r1.exists():
                continue

            # Per-read counts: [allele][edit_status] → int
            counts = defaultdict(lambda: defaultdict(int))
            n_read = n_processed = 0

            with utils.open_gz(r1) as f1, utils.open_gz(r2) as f2:
                for (hd1, s1, _, _), (_, s2, _, _) in zip(
                        utils.fastq_iter(f1), utils.fastq_iter(f2)):
                    n_read += 1
                    if args.max_reads > 0 and n_processed >= args.max_reads:
                        break

                    result = analyze_read(
                        s1, s2, anchor, allele_snps,
                        edit_offset, wt_base, desired_base,
                        min_snp_score_margin=args.min_allele_margin,
                        min_informative_snps=args.min_informative_snps,
                        edit_window_start=edit_win_start,
                        edit_window_end=edit_win_end,
                    )

                    al  = result["assigned_allele"]
                    es  = result["edit_status"]
                    counts[al][es] += 1

                    all_read_rows.append({
                        "target":                   t_name,
                        "sample_id":                sid,
                        "read_id":                  hd1.lstrip("@").split()[0],
                        "assigned_allele":          al,
                        "allele_score_chr01":       result["allele_score_chr01"],
                        "allele_score_chr11":       result["allele_score_chr11"],
                        "informative_snp_count":    result["informative_snp_count"],
                        "snp_base_observed":        result["snp_base_observed"],
                        "ambiguous_reason":         result["ambiguous_reason"],
                        "edit_position_base_obs":   result["edit_position_base_obs"],
                        "edit_status":              es,
                        "has_precise_desired_edit": result["has_precise_desired_edit"],
                    })
                    n_processed += 1

            print(f"[11_allele]   {t_name}/{sid}: {n_read} reads total, "
                  f"{n_processed} processed")

            # ── Build per-sample × per-allele summary ────────────────────
            m = sample_to_meta.get(sid, {})
            total_processed = n_processed

            # All assigned reads (chr01 + chr11), all categories
            for allele_name in allele_names + ["ambiguous", "low_information"]:
                ac        = counts[allele_name]
                allele_n  = sum(ac.values())
                wt_n      = ac["wt"]
                prec_n    = ac["precise_desired"]
                non_n     = ac["non_desired"]
                amb_edit  = ac["ambiguous_edit"]

                # Denominator 1: allele-assigned reads for this allele
                denom_allele     = allele_n
                pct_of_allele    = (prec_n / denom_allele * 100
                                    if denom_allele > 0 else 0.0)

                # Denominator 2: all informative allele-assigned reads
                total_assigned   = sum(
                    sum(counts[a].values()) for a in allele_names
                )
                pct_of_total_assigned = (prec_n / total_assigned * 100
                                         if total_assigned > 0 else 0.0)

                all_sample_rows.append({
                    "target":                     t_name,
                    "sample_id":                  sid,
                    "editor":                     m.get("editor", ""),
                    "dpi":                        m.get("dpi", ""),
                    "allele":                     allele_name,
                    "allele_assigned_reads":      allele_n,
                    "wt_reads":                   wt_n,
                    "precise_desired_reads":      prec_n,
                    "non_desired_reads":          non_n,
                    "ambiguous_edit_reads":       amb_edit,
                    "precise_desired_pct_of_allele_assigned":
                                                  f"{pct_of_allele:.4f}",
                    "precise_desired_pct_of_total_informative":
                                                  f"{pct_of_total_assigned:.4f}",
                    "VALIDATION_NOTE":            "NOT_VALIDATED_exploratory_only",
                })

            # QC row
            total_chr01   = sum(counts["chr01"].values())
            total_chr11   = sum(counts["chr11"].values())
            total_ambig   = sum(counts["ambiguous"].values())
            total_lowinfo = sum(counts["low_information"].values())
            assign_rate   = ((total_chr01 + total_chr11) /
                             total_processed * 100
                             if total_processed > 0 else 0.0)

            all_qc_rows.append({
                "target":                 t_name,
                "sample_id":              sid,
                "total_reads_processed":  total_processed,
                "chr01_reads":            total_chr01,
                "chr11_reads":            total_chr11,
                "ambiguous_reads":        total_ambig,
                "low_information_reads":  total_lowinfo,
                "allele_assignment_rate_pct": f"{assign_rate:.2f}",
                "VALIDATION_NOTE":        "NOT_VALIDATED_exploratory_only",
            })

    # ── Write primary outputs ────────────────────────────────────────────
    utils.write_csv(all_sample_rows,
                    out_dir / "allele_aware_precise_edit_summary.csv")
    utils.write_tsv(all_read_rows,
                    out_dir / "read_level_allele_edit_calls.tsv")
    utils.write_csv(all_qc_rows,
                    out_dir / "allele_assignment_qc.csv")

    print(f"[11_allele] Wrote {len(all_sample_rows)} summary rows to "
          f"allele_aware_precise_edit_summary.csv")
    print(f"[11_allele] Wrote {len(all_read_rows)} read-level rows to "
          f"read_level_allele_edit_calls.tsv")

    # ── Cross-check with validated motif_count_summary.tsv ──────────────
    motif_path = Path(args.input) / "04_edit" / "motif_count_summary.tsv"
    if motif_path.exists() and all_sample_rows:
        motif_rows = utils.read_tsv(motif_path)
        motif_lut  = {(r["target"], r["sample"]): r for r in motif_rows}

        compare_rows = []
        for sr in all_sample_rows:
            al = sr["allele"]
            if al not in ("chr01", "chr11"):
                continue

            # Map to corresponding allele motif-count target
            motif_target = f"ALSP_{al}"  # ALSP_chr01 or ALSP_chr11
            mk = (motif_target, sr["sample_id"])
            mr = motif_lut.get(mk)

            if mr:
                compare_rows.append({
                    "target":                  sr["target"],
                    "allele":                  al,
                    "sample_id":               sr["sample_id"],
                    "editor":                  sr["editor"],
                    "dpi":                     sr["dpi"],
                    # Allele-aware columns
                    "aa_allele_assigned_reads": sr["allele_assigned_reads"],
                    "aa_precise_desired_reads": sr["precise_desired_reads"],
                    "aa_precise_desired_pct":  sr["precise_desired_pct_of_allele_assigned"],
                    # Motif-count columns (validated)
                    "motif_total_pairs":       mr["total_read_pairs"],
                    "motif_desired_hits":      mr["desired_motif_hits"],
                    "motif_desired_pct":       mr["desired_percent_among_motif_hits"],
                    # Agreement
                    "VALIDATION_NOTE":         "aa_column_NOT_VALIDATED; motif_column_VALIDATED",
                })

        utils.write_csv(compare_rows,
                        out_dir / "allele_aware_vs_motif_comparison.csv")
        print(f"[11_allele] Wrote {len(compare_rows)} rows to "
              f"allele_aware_vs_motif_comparison.csv")
    else:
        # Ensure file exists even if empty
        utils.write_csv([{"note": "motif_count_summary.tsv not found or no allele-aware results"}],
                        out_dir / "allele_aware_vs_motif_comparison.csv")

    print(f"\n[11_allele] Outputs in {out_dir}")
    print(f"[11_allele] *** All outputs are exploratory and NOT validated ***")
    print(f"[11_allele] *** Validate against known controls before reporting ***")


if __name__ == "__main__":
    main()
