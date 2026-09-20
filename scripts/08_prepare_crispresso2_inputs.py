#!/usr/bin/env python3
"""
08_prepare_crispresso2_inputs.py  --  Prepare CRISPResso2 batch input files.

*** OPTIONAL MODULE — NOT PART OF THE VALIDATED MOTIF-COUNT PIPELINE ***

Creates CRISPResso2-compatible per-sample FASTQs and a CRISPRessoBatch CSV
from demultiplexed HiTOM reads. CRISPResso2 must be installed separately;
this script only prepares the inputs.

LIMITATIONS:
  - amplicon_seq in the batch CSV requires a full reference amplicon sequence.
  - If reference_amplicon_sequence = null in targets.json, this field is left
    as a placeholder that must be filled in manually.
  - expected_hdr_amplicon_seq similarly requires the full edited amplicon.
  - Do not run CRISPResso2 on these files until amplicon sequences are verified.

Outputs (under <output>/08_crispresso2_inputs/):
  fastqs/           symlinks or copies of demux_clean per-sample FASTQs
  CRISPRessoBatch_input.csv
  README_crispresso2_inputs.md
"""

import argparse
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils

PIPELINE_ROOT    = Path(__file__).resolve().parent.parent
DEFAULT_METADATA = PIPELINE_ROOT / "config" / "sample_metadata.csv"
DEFAULT_CONFIG   = PIPELINE_ROOT / "config" / "targets.json"

PLACEHOLDER = "FILL_IN_AMPLICON_SEQUENCE"


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input",    required=True,
                   help="Root output dir (expects <input>/01_demux/demux_clean/)")
    p.add_argument("--output",   required=True,
                   help="Root output dir (outputs to <output>/08_crispresso2_inputs/)")
    p.add_argument("--metadata", default=str(DEFAULT_METADATA))
    p.add_argument("--config",   default=str(DEFAULT_CONFIG))
    p.add_argument("--use-symlinks", action="store_true", dest="symlinks",
                   help="Symlink FASTQs instead of copying (saves disk space)")
    return p.parse_args()


def main():
    args      = parse_args()
    clean_dir = Path(args.input) / "01_demux" / "demux_clean"
    out_dir   = Path(args.output) / "08_crispresso2_inputs"
    fastq_dir = out_dir / "fastqs"

    for d in [out_dir, fastq_dir]:
        d.mkdir(parents=True, exist_ok=True)

    if not clean_dir.exists():
        print(f"[08_crispresso2] ERROR: demux_clean not found: {clean_dir}",
              file=sys.stderr)
        sys.exit(1)

    meta_rows       = utils.read_csv_file(args.metadata)
    targets_cfg     = utils.load_targets(args.config)["targets"]
    target_by_name  = {t["name"]: t for t in targets_cfg}
    # Map target_gene → list of target entries
    target_by_gene  = {}
    for t in targets_cfg:
        for app in t["applicable_to"]:
            target_by_gene.setdefault(app.upper(), []).append(t)

    batch_rows = []
    missing_amplicons = []

    for row in meta_rows:
        sid        = row["sample_id"]
        tgene      = row["target"].upper()
        r1_src     = clean_dir / f"{sid}_R1.fastq.gz"
        r2_src     = clean_dir / f"{sid}_R2.fastq.gz"

        if not r1_src.exists():
            print(f"[08_crispresso2] SKIP {sid}: R1 not found in demux_clean")
            continue

        # Symlink or copy
        r1_dst = fastq_dir / r1_src.name
        r2_dst = fastq_dir / r2_src.name
        for src, dst in [(r1_src, r1_dst), (r2_src, r2_dst)]:
            if not dst.exists():
                if args.symlinks:
                    os.symlink(src.resolve(), dst)
                else:
                    import shutil
                    shutil.copy2(src, dst)

        # One CRISPResso2 row per target applicable to this sample
        for tgt in target_by_gene.get(tgene, []):
            amp_seq     = tgt.get("reference_amplicon_sequence") or PLACEHOLDER
            hdr_seq     = tgt.get("expected_edited_amplicon_sequence") or ""
            guide_seq   = tgt.get("guide_sequence") or ""

            if amp_seq == PLACEHOLDER:
                missing_amplicons.append((sid, tgt["name"]))

            name_tag = f"{sid}__{tgt['name']}"
            batch_rows.append({
                "name":                      name_tag,
                "fastq_r1":                  str(r1_dst.resolve()),
                "fastq_r2":                  str(r2_dst.resolve()),
                "amplicon_seq":              amp_seq,
                "expected_hdr_amplicon_seq": hdr_seq,
                "guide_seq":                 guide_seq,
                "target_name":               tgt["name"],
                "sample_id":                 sid,
                "editor":                    row.get("editor", ""),
                "dpi":                       row.get("dpi", ""),
                "_NOTE":                     (
                    "amplicon_seq=PLACEHOLDER — fill in before running CRISPResso2"
                    if amp_seq == PLACEHOLDER else ""),
            })

    # Write batch CSV (CRISPResso2-compatible columns first)
    c2_cols = ["name", "fastq_r1", "fastq_r2",
               "amplicon_seq", "expected_hdr_amplicon_seq", "guide_seq"]
    extra   = ["target_name", "sample_id", "editor", "dpi", "_NOTE"]

    with open(out_dir / "CRISPRessoBatch_input.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=c2_cols + extra)
        w.writeheader()
        w.writerows(batch_rows)

    # README
    readme_lines = [
        "# CRISPResso2 Batch Input Files",
        "",
        "Generated by `08_prepare_crispresso2_inputs.py`.",
        "",
        "## Before running CRISPResso2",
        "",
    ]
    if missing_amplicons:
        readme_lines += [
            f"**{len(missing_amplicons)} amplicon sequences are missing** "
            f"(marked as FILL_IN_AMPLICON_SEQUENCE).",
            "",
            "Edit `CRISPRessoBatch_input.csv` and replace "
            "`FILL_IN_AMPLICON_SEQUENCE` with the correct full amplicon sequence",
            "for each target before running CRISPResso2.",
            "",
            "Samples/targets needing amplicon sequences:",
        ] + [f"  - {sid} / {tname}" for sid, tname in missing_amplicons] + [""]
    else:
        readme_lines += ["All amplicon sequences are filled in. Ready to run.", ""]

    readme_lines += [
        "## Run command",
        "",
        "```bash",
        "CRISPRessoBatch \\",
        "    --batch_settings CRISPRessoBatch_input.csv \\",
        "    --output_folder ../09_crispresso2_results \\",
        "    --min_frequency_alleles_around_cut_to_plot 0.05 \\",
        "    -n 4",
        "```",
        "",
        "See `pipeline/docs/crispresso2_local_usage.md` for installation "
        "and interpretation.",
    ]

    (out_dir / "README_crispresso2_inputs.md").write_text("\n".join(readme_lines))

    print(f"[08_crispresso2] Wrote {len(batch_rows)} batch rows "
          f"to CRISPRessoBatch_input.csv")
    if missing_amplicons:
        print(f"[08_crispresso2] WARNING: {len(missing_amplicons)} entries have "
              f"FILL_IN_AMPLICON_SEQUENCE — fill in before running CRISPResso2")
    print(f"[08_crispresso2] Outputs in {out_dir}")


if __name__ == "__main__":
    main()
