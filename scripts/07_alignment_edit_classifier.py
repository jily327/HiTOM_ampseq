#!/usr/bin/env python3
"""
07_alignment_edit_classifier.py  --  Alignment-based edit classification.

*** OPTIONAL MODULE — NOT PART OF THE VALIDATED MOTIF-COUNT PIPELINE ***
*** ALL OUTPUTS ARE NEWLY IMPLEMENTED AND REQUIRE VALIDATION          ***

Uses local pairwise alignment (Smith-Waterman, pure Python) to classify
reads from the HiTOM demuxed FASTQs. This extends the motif-count pipeline
by detecting indels, imprecise prime edits, and other substitutions in the
motif window.

IMPORTANT LIMITATIONS:
  - Because full amplicon reference sequences are not available for this
    dataset (reference_amplicon_sequence = null in targets.json), alignment
    is performed against the short WT/desired motif sequences (~15-29 bp)
    rather than the full amplicon.
  - This means indels and substitutions OUTSIDE the motif window are NOT
    detected. Results represent the motif window only.
  - Reads already captured by the validated motif-count pipeline (exact
    WT or desired match) are still re-classified here for completeness,
    but the validated motif_count_summary.tsv is NEVER modified.
  - These outputs have NOT been validated against known controls. Treat
    as exploratory until confirmed.

Classification categories:
  wt               Exact WT motif match (consistent with motif-count pipeline)
  precise_desired  Exact desired motif match, no other variants in motif window
  imprecise_PE     Desired-motif-like but with additional mismatches
  indel            Indel detected in motif window, no exact motif match
  desired_plus_indel  Desired motif match AND indel detected in alignment
  other_substitution  Aligned to WT with substitution(s) but not the desired change
  other_allele     Exact match to a different allele's motif of the same
                   gene (e.g. a chr11 read scored against the chr01 target).
                   Reported separately so that it is not mistaken for an
                   imprecise edit or a failed alignment.
  no_motif_align   No alignment to either motif (score below threshold)
  short_read       Read too short for reliable alignment (< --min-read-len)

READ SAMPLING:
  --max-reads limits how many read pairs per sample/target are aligned
  (Smith-Waterman in pure Python is slow).  The subset is drawn at random
  from the whole FASTQ by default (--sampling random, reservoir sampling with
  --seed), not taken from the head of the file, so it is not biased by read
  order on the flowcell.  --sampling head restores the previous first-N
  behaviour.  Every summary row records the full read count of the sample,
  how many reads were classified, and the strategy and seed used.

Outputs (under <output>/07_alignment/):
  alignment_read_classifications.tsv
  alignment_sample_summary.csv
  alignment_indel_summary.csv
  alignment_window_mutation_spectrum.csv
  alignment_vs_motif_comparison.csv
"""

import argparse
import csv
import random
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils

PIPELINE_ROOT    = Path(__file__).resolve().parent.parent
DEFAULT_METADATA = PIPELINE_ROOT / "config" / "sample_metadata.csv"
DEFAULT_CONFIG   = PIPELINE_ROOT / "config" / "targets.json"

# Alignment parameters
MATCH   =  2
MISMATCH = -1
GAP      = -2
MIN_IDENTITY_DEFAULT = 0.70   # minimum alignment identity to call a classification
MAX_READS_DEFAULT    = 3000   # reads per sample (per target) for speed; 0 = all
SAMPLING_DEFAULT     = "random"  # how the --max-reads subset is drawn
SEED_DEFAULT         = 0         # RNG seed, so a run is reproducible


# ---------------------------------------------------------------------------
# Pure-Python Smith-Waterman local alignment
# ---------------------------------------------------------------------------

def sw_align(query: str, ref: str,
             match: int = MATCH,
             mismatch: int = MISMATCH,
             gap: int = GAP):
    """
    Smith-Waterman local alignment.

    Returns: (score, aligned_query, aligned_ref, ref_coverage_identity)

    ref_coverage_identity = matching_positions / len(ref)
      This tells us how much of the reference motif is matched by the query.
    """
    q = query.upper()
    r = ref.upper()
    n, m = len(q), len(r)

    # DP matrix (n+1 × m+1), row-major list-of-lists
    H = [[0] * (m + 1) for _ in range(n + 1)]
    best_score = 0
    best_i, best_j = 0, 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            sc = match if q[i - 1] == r[j - 1] else mismatch
            cell = max(0,
                       H[i - 1][j - 1] + sc,
                       H[i - 1][j]     + gap,
                       H[i][j - 1]     + gap)
            H[i][j] = cell
            if cell > best_score:
                best_score = cell
                best_i, best_j = i, j

    if best_score == 0:
        return 0, "", "", 0.0

    # Traceback
    aq, ar = [], []
    i, j = best_i, best_j
    while i > 0 and j > 0 and H[i][j] > 0:
        sc = match if q[i - 1] == r[j - 1] else mismatch
        if H[i][j] == H[i - 1][j - 1] + sc:
            aq.append(q[i - 1])
            ar.append(r[j - 1])
            i -= 1; j -= 1
        elif H[i][j] == H[i - 1][j] + gap:
            aq.append(q[i - 1])
            ar.append('-')
            i -= 1
        else:
            aq.append('-')
            ar.append(r[j - 1])
            j -= 1

    aq = ''.join(reversed(aq))
    ar = ''.join(reversed(ar))

    # Identity relative to reference length (how much of motif was covered)
    matches = sum(a == b for a, b in zip(aq, ar) if a != '-' and b != '-')
    ref_cov_id = matches / m if m > 0 else 0.0

    return best_score, aq, ar, ref_cov_id


def find_edit_positions(wt_motif: str, desired_motif: str):
    """
    Return list of (position, wt_base, desired_base) from WT/desired diff.
    Positions are 0-indexed within the motif.
    """
    return [
        (i, wt_motif[i].upper(), desired_motif[i].upper())
        for i in range(min(len(wt_motif), len(desired_motif)))
        if wt_motif[i].upper() != desired_motif[i].upper()
    ]


def sibling_allele_motifs(target: dict, all_targets: list):
    """
    Return [(allele_name, motif_seq), ...] for the OTHER alleles of the same
    gene, so that reads coming from the wrong homeolog can be labelled
    'other_allele' instead of being scored against this allele's motifs.

    A sibling is a target that covers the same gene(s) (identical
    applicable_to set) and carries a different, non-null allele name.
    Returns an empty list for allele-agnostic targets.
    """
    if not target.get("allele"):
        return []
    app = {a.upper() for a in target["applicable_to"]}
    out = []
    for other in all_targets:
        if other is target or not other.get("allele"):
            continue
        if other["allele"] == target["allele"]:
            continue
        if {a.upper() for a in other["applicable_to"]} != app:
            continue
        out.append((other["allele"], other["wt_motif"].upper()))
        out.append((other["allele"], other["desired_motif"].upper()))
    return out


def detect_indels_in_alignment(aligned_q: str, aligned_r: str):
    """
    Scan a pairwise alignment for insertion and deletion events.
    Returns: list of (type, position_in_ref, size)
    type: 'del' (gap in query) or 'ins' (gap in ref)
    """
    events = []
    ref_pos = 0
    in_del = in_ins = False
    del_start = ins_start = 0
    del_size = ins_size = 0

    for i, (q_c, r_c) in enumerate(zip(aligned_q, aligned_r)):
        if r_c == '-':   # insertion in query
            if not in_ins:
                in_ins = True
                ins_start = ref_pos
                ins_size = 0
            ins_size += 1
            if in_del:
                events.append(('del', del_start, del_size))
                in_del = False
        elif q_c == '-': # deletion in query (gap)
            if not in_del:
                in_del = True
                del_start = ref_pos
                del_size = 0
            del_size += 1
            ref_pos += 1
            if in_ins:
                events.append(('ins', ins_start, ins_size))
                in_ins = False
        else:
            if in_del:
                events.append(('del', del_start, del_size))
                in_del = False
            if in_ins:
                events.append(('ins', ins_start, ins_size))
                in_ins = False
            ref_pos += 1

    if in_del:
        events.append(('del', del_start, del_size))
    if in_ins:
        events.append(('ins', ins_start, ins_size))

    return events


def classify_read(s1: str, s2: str,
                  wt_motif: str, desired_motif: str,
                  edit_positions: list,
                  min_identity: float = MIN_IDENTITY_DEFAULT,
                  min_read_len: int = 30,
                  other_allele_motifs=()):
    """
    Classify one read pair (s1=R1, s2=R2) against WT and desired motifs.

    other_allele_motifs: [(allele_name, motif_seq), ...] from the other
    homeolog(s) of the same gene.  A read that matches one of them exactly,
    and neither of this target's own motifs, is reported as 'other_allele'.

    Returns: dict with classification fields.
    """
    combined = s1.upper() + "N" + s2.upper() + "N" + utils.reverse_complement(s1) + "N" + utils.reverse_complement(s2)
    read_len = max(len(s1), len(s2))

    if read_len < min_read_len:
        return {"classification": "short_read",
                "best_reference": "none",
                "alignment_score": 0,
                "alignment_identity": 0.0,
                "has_desired_edit": False,
                "has_indel": False,
                "has_mismatch": False,
                "indel_size": 0,
                "indel_near_edit_window": False}

    # -- Exact motif check (same logic as validated 04_edit_classification.py) --
    exact_wt      = utils.motif_hit(combined, wt_motif,      max_mismatches=0)
    exact_desired = utils.motif_hit(combined, desired_motif, max_mismatches=0)

    if exact_wt and not exact_desired:
        return {"classification": "wt",
                "best_reference": "wt_motif",
                "alignment_score": len(wt_motif) * MATCH,
                "alignment_identity": 1.0,
                "has_desired_edit": False,
                "has_indel": False,
                "has_mismatch": False,
                "indel_size": 0,
                "indel_near_edit_window": False}

    if exact_desired and not exact_wt:
        # Verify no indels in the surrounding sequence via alignment
        score, aq, ar, ident = sw_align(combined, desired_motif)
        indels = detect_indels_in_alignment(aq, ar) if aq else []
        has_indel = len(indels) > 0
        total_indel_size = sum(e[2] for e in indels)
        return {"classification": "desired_plus_indel" if has_indel else "precise_desired",
                "best_reference": "desired_motif",
                "alignment_score": score,
                "alignment_identity": ident,
                "has_desired_edit": True,
                "has_indel": has_indel,
                "has_mismatch": False,
                "indel_size": total_indel_size,
                "indel_near_edit_window": has_indel}

    if not exact_wt and not exact_desired:
        # The read may simply belong to the other homeolog.  Without this
        # check such reads are aligned against this allele's motifs and land
        # in imprecise_PE / other_substitution / no_motif_align, which reads
        # as an editing outcome when it is only allele cross-talk.
        for allele_name, other_motif in other_allele_motifs:
            if utils.motif_hit(combined, other_motif, max_mismatches=0):
                return {"classification": "other_allele",
                        "best_reference": f"{allele_name}_motif",
                        "alignment_score": len(other_motif) * MATCH,
                        "alignment_identity": 1.0,
                        "has_desired_edit": False,
                        "has_indel": False,
                        "has_mismatch": False,
                        "indel_size": 0,
                        "indel_near_edit_window": False,
                        "other_allele_match": allele_name}

    if exact_wt and exact_desired:
        return {"classification": "precise_desired",
                "best_reference": "both_motifs",
                "alignment_score": len(desired_motif) * MATCH,
                "alignment_identity": 1.0,
                "has_desired_edit": True,
                "has_indel": False,
                "has_mismatch": False,
                "indel_size": 0,
                "indel_near_edit_window": False}

    # -- No exact match: try alignment to both references --
    score_wt,  aq_wt,  ar_wt,  id_wt  = sw_align(combined, wt_motif)
    score_des, aq_des, ar_des, id_des  = sw_align(combined, desired_motif)

    use_wt    = id_wt  >= min_identity
    use_des   = id_des >= min_identity

    if not use_wt and not use_des:
        return {"classification": "no_motif_align",
                "best_reference": "none",
                "alignment_score": max(score_wt, score_des),
                "alignment_identity": max(id_wt, id_des),
                "has_desired_edit": False,
                "has_indel": False,
                "has_mismatch": False,
                "indel_size": 0,
                "indel_near_edit_window": False}

    # Use the better alignment
    if score_des >= score_wt and use_des:
        best_ref  = "desired_motif"
        aq, ar    = aq_des, ar_des
        best_id   = id_des
        best_sc   = score_des
    else:
        best_ref  = "wt_motif"
        aq, ar    = aq_wt, ar_wt
        best_id   = id_wt
        best_sc   = score_wt

    indels       = detect_indels_in_alignment(aq, ar)
    has_indel    = len(indels) > 0
    total_size   = sum(e[2] for e in indels)

    # Check if any indel overlaps with the edit window
    edit_pos_set = {p for p, _, _ in edit_positions}
    near_edit    = any(
        any(ep in range(ev_pos, ev_pos + ev_size + 1) for ep in edit_pos_set)
        for ev_type, ev_pos, ev_size in indels
    )

    # Check for imprecise PE: desired-like changes near edit site but not exact
    mismatches_at_edit = []
    if aq and ar and best_ref == "wt_motif":
        ref_offset = 0
        for qc, rc in zip(aq, ar):
            if rc == '-':
                continue
            if qc != '-' and qc != rc:
                mismatches_at_edit.append((ref_offset, rc, qc))
            if rc != '-':
                ref_offset += 1

    desired_changes_at_edit = set((p, wt, des) for p, wt, des in edit_positions)
    observed_at_edit = set((pos, rc.upper(), qc.upper())
                           for pos, rc, qc in mismatches_at_edit
                           if (pos, rc.upper(), qc.upper()) in desired_changes_at_edit)

    is_imprecise = (len(mismatches_at_edit) > 0 and
                    len(observed_at_edit) < len(desired_changes_at_edit) and
                    not has_indel)
    is_other_sub = (len(mismatches_at_edit) > 0 and
                    len(observed_at_edit) == 0 and
                    not has_indel)

    if has_indel and best_ref == "desired_motif":
        classification = "desired_plus_indel"
    elif has_indel:
        classification = "indel"
    elif is_imprecise:
        classification = "imprecise_PE"
    elif is_other_sub:
        classification = "other_substitution"
    else:
        classification = "no_motif_align"

    return {"classification": classification,
            "best_reference": best_ref,
            "alignment_score": best_sc,
            "alignment_identity": round(best_id, 4),
            "has_desired_edit": best_ref == "desired_motif",
            "has_indel": has_indel,
            "has_mismatch": len(mismatches_at_edit) > 0,
            "indel_size": total_size,
            "indel_near_edit_window": near_edit}


# ---------------------------------------------------------------------------
# Read selection
# ---------------------------------------------------------------------------

def stream_read_pairs(r1_path, r2_path, label=""):
    """Yield (read_id, s1, s2) for every pair in a demuxed FASTQ pair."""
    with utils.open_gz(r1_path) as f1, utils.open_gz(r2_path) as f2:
        for (hd1, s1, _, _), (_, s2, _, _) in utils.fastq_pair_iter(f1, f2, label):
            yield hd1.lstrip("@").split()[0], s1, s2


def select_read_pairs(r1_path, r2_path, max_reads, strategy, rng, label=""):
    """
    Return (selected_pairs, total_pairs_in_sample).

    strategy 'random' draws an unbiased subset of size max_reads from the
    whole file by reservoir sampling; 'head' keeps the first max_reads pairs.
    Either way the file is read to the end so that the reported total read
    count is the real one, not just the number that happened to be scored.
    Memory is bounded by max_reads, not by the size of the FASTQ.
    """
    selected = []
    total = 0
    for rec in stream_read_pairs(r1_path, r2_path, label):
        total += 1
        if len(selected) < max_reads:
            selected.append(rec)
        elif strategy == "random":
            j = rng.randrange(total)
            if j < max_reads:
                selected[j] = rec
    return selected, total


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
                   help="Root output dir (outputs go to <output>/07_alignment/)")
    p.add_argument("--metadata", default=str(DEFAULT_METADATA))
    p.add_argument("--config",   default=str(DEFAULT_CONFIG))
    p.add_argument("--target",   default=None,
                   help="Restrict to one target name (e.g. ALSW). Default: all targets.")
    p.add_argument("--min-map-identity", type=float, default=MIN_IDENTITY_DEFAULT,
                   dest="min_identity",
                   help=f"Minimum ref-coverage identity to call alignment "
                        f"[default: {MIN_IDENTITY_DEFAULT}]")
    p.add_argument("--min-read-len", type=int, default=30, dest="min_read_len",
                   help="Minimum R1 length to attempt alignment [default: 30]")
    p.add_argument("--max-reads", type=int, default=MAX_READS_DEFAULT, dest="max_reads",
                   help="Max read pairs per sample per target (0=all) "
                        f"[default: {MAX_READS_DEFAULT}]")
    p.add_argument("--sampling", choices=["random", "head"], default=SAMPLING_DEFAULT,
                   help="How the --max-reads subset is drawn: 'random' "
                        "(reservoir sampling over the whole file, unbiased) or "
                        "'head' (first N reads, biased by read order) "
                        f"[default: {SAMPLING_DEFAULT}]")
    p.add_argument("--seed", type=int, default=SEED_DEFAULT,
                   help=f"RNG seed for --sampling random [default: {SEED_DEFAULT}]")
    p.add_argument("--window-start", type=int, default=None, dest="window_start",
                   help="Override edit window start (0-indexed, default: from targets.json)")
    p.add_argument("--window-end",   type=int, default=None, dest="window_end",
                   help="Override edit window end (0-indexed, default: from targets.json)")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args      = parse_args()
    clean_dir = Path(args.input) / "01_demux" / "demux_clean"
    out_dir   = Path(args.output) / "07_alignment"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not clean_dir.exists():
        print(f"[07_align] ERROR: demux_clean not found: {clean_dir}", file=sys.stderr)
        print(f"[07_align]   Run 01_demux_and_qc.py first.", file=sys.stderr)
        sys.exit(1)

    meta_rows = utils.read_csv_file(args.metadata)
    sample_to_meta = {r["sample_id"]: r for r in meta_rows}

    targets_cfg = utils.load_targets(args.config)["targets"]
    if args.target:
        targets_cfg = [t for t in targets_cfg if t["name"] == args.target]
        if not targets_cfg:
            print(f"[07_align] ERROR: target '{args.target}' not in targets.json",
                  file=sys.stderr)
            sys.exit(1)

    print(f"[07_align] *** NOT VALIDATED — exploratory analysis only ***")
    print(f"[07_align] NOTE: alignment uses short WT/desired motifs as reference.")
    print(f"[07_align]       Full amplicon alignment requires reference_amplicon_sequence.")
    print(f"[07_align] Targets: {[t['name'] for t in targets_cfg]}")
    if args.max_reads == 0:
        sampling_label = "all"
        print(f"[07_align] Max reads per sample/target: all")
    else:
        sampling_label = args.sampling
        print(f"[07_align] Max reads per sample/target: {args.max_reads} "
              f"({args.sampling} sampling, seed={args.seed})")
        if args.sampling == "head":
            print(f"[07_align] WARNING: --sampling head takes the first "
                  f"{args.max_reads} reads of each file. That subset is "
                  f"ordered by flowcell position, so it is not a random "
                  f"sample of the library.")
    rng = random.Random(args.seed)

    # Sibling allele motifs are resolved from the FULL target list, so that
    # --target still knows about the other homeolog.
    all_targets = utils.load_targets(args.config)["targets"]

    read_class_rows       = []
    sample_summary        = defaultdict(lambda: defaultdict(int))
    sample_reads          = {}   # (target, sample) → read pairs in the FASTQ
    indel_rows            = []
    mutation_spectrum     = []

    for target in targets_cfg:
        tname         = target["name"]
        wt_motif      = target["wt_motif"].upper()
        desired_motif = target["desired_motif"].upper()
        applicable    = [a.upper() for a in target["applicable_to"]]
        edit_pos      = find_edit_positions(wt_motif, desired_motif)
        other_motifs  = sibling_allele_motifs(target, all_targets)
        if other_motifs:
            print(f"[07_align] {tname}: reads matching "
                  f"{sorted({a for a, _ in other_motifs})} motifs will be "
                  f"reported as 'other_allele'")

        # Warn if using motif as short reference
        if target.get("reference_amplicon_sequence") is None:
            print(f"[07_align] {tname}: no full amplicon ref → using "
                  f"{len(wt_motif)}bp motif as reference (motif-window only)")

        # Iterate samples in metadata order
        for row in meta_rows:
            sid    = row["sample_id"]
            if row["target"].upper() not in applicable:
                continue
            r1 = clean_dir / f"{sid}_R1.fastq.gz"
            r2 = clean_dir / f"{sid}_R2.fastq.gz"
            if not r1.exists():
                continue

            label = f"{tname} / {sid}"
            if args.max_reads > 0:
                pairs, n_read = select_read_pairs(
                    r1, r2, args.max_reads, args.sampling, rng, label)
            else:
                pairs, n_read = stream_read_pairs(r1, r2, label), None

            n_short = n_processed = 0
            for read_id, s1, s2 in pairs:
                result = classify_read(
                    s1, s2, wt_motif, desired_motif,
                    edit_pos, args.min_identity, args.min_read_len,
                    other_allele_motifs=other_motifs,
                )
                clss   = result["classification"]
                sample_summary[(tname, sid)][clss] += 1

                if clss == "short_read":
                    n_short += 1
                n_processed += 1

                read_class_rows.append({
                    "target":                   tname,
                    "sample_id":                sid,
                    "read_id":                  read_id,
                    "classification":           clss,
                    "best_reference":           result["best_reference"],
                    "alignment_score":          result["alignment_score"],
                    "alignment_identity":       result["alignment_identity"],
                    "has_desired_edit":         result["has_desired_edit"],
                    "has_indel":                result["has_indel"],
                    "indel_size":               result["indel_size"],
                    "indel_near_edit_window":   result["indel_near_edit_window"],
                    "other_allele_match":       result.get("other_allele_match", ""),
                })

                if result["has_indel"] and clss != "short_read":
                    indel_rows.append({
                        "target":        tname,
                        "sample_id":     sid,
                        "read_id":       read_id,
                        "classification": clss,
                        "indel_size":    result["indel_size"],
                        "near_edit":     result["indel_near_edit_window"],
                    })

            if n_read is None:          # streamed every read
                n_read = n_processed
            sample_reads[(tname, sid)] = n_read

            print(f"[07_align]   {tname} / {sid}: {n_read} reads in sample, "
                  f"{n_processed} classified ({n_short} short)")

    # -----------------------------------------------------------------------
    # Build sample summary
    # -----------------------------------------------------------------------
    CATS = ["wt", "precise_desired", "imprecise_PE", "indel",
            "desired_plus_indel", "other_substitution", "other_allele",
            "no_motif_align", "short_read"]

    sample_summary_rows = []
    for (tname, sid), counts in sample_summary.items():
        m         = sample_to_meta.get(sid, {})
        total     = sum(counts.values())
        informative = counts["wt"] + counts["precise_desired"]
        prec_pct  = (counts["precise_desired"] / informative * 100
                     if informative > 0 else 0.0)
        in_sample = sample_reads.get((tname, sid), total)
        frac      = (total / in_sample * 100) if in_sample else 0.0
        row = {
            "target":   tname,
            "sample_id": sid,
            "editor":   m.get("editor", ""),
            "dpi":      m.get("dpi", ""),
            "total_reads_processed": total,
            # Read accounting: how many reads the sample actually has, how many
            # of them were classified here, and how the subset was drawn.  The
            # motif-count pipeline (04) always uses every read, so a comparison
            # between the two is only meaningful alongside these columns.
            "total_read_pairs_in_sample": in_sample,
            "reads_classified_pct": f"{frac:.2f}",
            "sampling_strategy": ("all" if args.max_reads == 0
                                  else args.sampling),
            "sampling_seed": ("" if args.max_reads == 0 or args.sampling == "head"
                              else args.seed),
            "alignment_precise_desired_pct": f"{prec_pct:.4f}",
            "VALIDATION_NOTE": "NOT_VALIDATED_exploratory_only",
        }
        for cat in CATS:
            row[cat] = counts.get(cat, 0)
        sample_summary_rows.append(row)

    # -----------------------------------------------------------------------
    # Write outputs
    # -----------------------------------------------------------------------
    utils.write_tsv(read_class_rows,
                    out_dir / "alignment_read_classifications.tsv")
    utils.write_csv(sample_summary_rows,
                    out_dir / "alignment_sample_summary.csv")
    utils.write_csv(indel_rows or [{"note": "no indels detected"}],
                    out_dir / "alignment_indel_summary.csv")

    # Alignment vs motif comparison (join with motif_count_summary.tsv)
    motif_path = Path(args.input) / "04_edit" / "motif_count_summary.tsv"
    if motif_path.exists():
        motif_rows = utils.read_tsv(motif_path)
        motif_lut  = {(r["target"], r["sample"]): r for r in motif_rows}
        compare    = []
        for srow in sample_summary_rows:
            mk = (srow["target"], srow["sample_id"])
            if mk in motif_lut:
                mr = motif_lut[mk]
                compare.append({
                    "target":            srow["target"],
                    "sample_id":         srow["sample_id"],
                    "motif_desired_pct": mr["desired_percent_among_motif_hits"],
                    "alignment_desired_pct": srow["alignment_precise_desired_pct"],
                    "motif_informative": mr["informative_reads"],
                    "motif_total_read_pairs": mr["total_read_pairs"],
                    "alignment_total_read_pairs": srow["total_read_pairs_in_sample"],
                    "alignment_reads_classified": srow["total_reads_processed"],
                    "alignment_sampling": srow["sampling_strategy"],
                    "alignment_wt":      srow["wt"],
                    "alignment_precise": srow["precise_desired"],
                    "alignment_indel":   srow["indel"],
                    "alignment_imprecise": srow["imprecise_PE"],
                    "alignment_other_allele": srow["other_allele"],
                    "VALIDATION_NOTE":   "alignment_column_NOT_VALIDATED",
                })
        utils.write_csv(compare,
                        out_dir / "alignment_vs_motif_comparison.csv")

    # Mutation spectrum placeholder
    utils.write_csv([{"note": "mutation_spectrum_requires_full_amplicon_reference"}],
                    out_dir / "alignment_window_mutation_spectrum.csv")

    print(f"\n[07_align] Outputs in {out_dir}")
    print(f"[07_align] *** All outputs are exploratory and NOT validated ***")
    print(f"[07_align] *** Do NOT use these values in publications without ")
    print(f"[07_align]     independent validation against known controls   ***")


if __name__ == "__main__":
    main()
