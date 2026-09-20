#!/usr/bin/env python3
"""
make_new_hitom_project.py  --  Generate a new HiTOM pipeline project.

Given a CSV listing sample IDs and barcode assignments, creates a complete
ready-to-run project directory with config files, an example run command,
and an empty expected output structure.

This lets you reuse the pipeline for future HiTOM runs without editing any
Python code. You only need to:
  1. Provide a sample sheet CSV.
  2. Fill in target motif sequences in the generated targets.json.
  3. Run the pipeline.

Input CSV must have at minimum these columns:
  sample_id, target, editor, dpi, forward_barcode, reverse_barcode

Optional columns: leaf, replicate, notes

Barcode columns accept either HiTOM names (F1, R-A) or raw 4-bp sequences.

Usage:
  python3 make_new_hitom_project.py \\
      --samples my_new_samples.csv \\
      --project-dir ~/projects/my_new_hitom_run \\
      --project-name "My New Experiment"
"""

import argparse
import csv
import json
import shutil
import sys
import textwrap
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PROTOCOL = PIPELINE_ROOT / "config" / "hitom_protocol_constants.json"

REQUIRED_COLS = {"sample_id", "target", "editor", "dpi",
                 "forward_barcode", "reverse_barcode"}

TARGETS_JSON_TEMPLATE = {
    "_instructions": (
        "Fill in the fields below for each target in your experiment. "
        "wt_motif and desired_motif are REQUIRED for motif-count analysis. "
        "reference_amplicon_sequence is OPTIONAL (for alignment-based analysis). "
        "See pipeline/config/targets.json from the ALSW project for a worked example."
    ),
    "targets": []
}

TARGET_ENTRY_TEMPLATE = {
    "name":                        "TARGET_NAME_HERE",
    "applicable_to":               ["TARGET_GENE_HERE"],
    "wt_motif":                    "FILL_IN_WT_MOTIF_SEQUENCE",
    "desired_motif":               "FILL_IN_DESIRED_MOTIF_SEQUENCE",
    "allele":                      None,
    "category":                    "edit_classification",
    "note":                        "Fill in notes for this target.",
    "edit_positions_in_motif":     [],
    "reference_amplicon_sequence": None,
    "_reference_amplicon_note":    "Optional. Provide the full amplicon sequence for alignment-based analysis.",
    "expected_edited_amplicon_sequence": None,
    "guide_sequence":              None,
    "analysis_notes":              ""
}


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--samples", required=True,
                   help="CSV file with sample sheet (sample_id, target, editor, dpi, "
                        "forward_barcode, reverse_barcode, [leaf, replicate])")
    p.add_argument("--project-dir", required=True, dest="project_dir",
                   help="Path for new project directory (will be created)")
    p.add_argument("--project-name", default=None, dest="project_name",
                   help="Human-readable project name [default: directory name]")
    p.add_argument("--protocol-constants", default=str(DEFAULT_PROTOCOL),
                   dest="protocol_constants",
                   help="Path to hitom_protocol_constants.json "
                        "[default: pipeline/config/hitom_protocol_constants.json]")
    p.add_argument("--copy-pipeline", default=str(PIPELINE_ROOT),
                   dest="copy_pipeline",
                   help="Pipeline source to copy into project [default: current pipeline]")
    p.add_argument("--raw-fastq-dir", default=None, dest="raw_fastq_dir",
                   help="Optional: path to raw FASTQ directory (written into example command)")
    return p.parse_args()


def load_protocol(path):
    with open(path) as f:
        return json.load(f)


def validate_sample_sheet(rows, fw_dict, rv_dict):
    """Check that all required columns exist and all barcodes are known."""
    if not rows:
        raise ValueError("Sample sheet is empty.")

    cols = set(rows[0].keys())
    missing = REQUIRED_COLS - cols
    if missing:
        raise ValueError(f"Sample sheet missing required columns: {missing}")

    errors = []
    sample_ids = []
    for row in rows:
        sid = row["sample_id"]
        if sid in sample_ids:
            errors.append(f"Duplicate sample_id: {sid}")
        sample_ids.append(sid)

        # Validate barcodes
        fb = row["forward_barcode"].strip().upper()
        rb = row["reverse_barcode"].strip().upper()

        if fb not in fw_dict and not all(c in "ACGT" for c in fb):
            errors.append(f"{sid}: unknown forward_barcode '{fb}'")
        if rb not in rv_dict and not all(c in "ACGT" for c in rb):
            errors.append(f"{sid}: unknown reverse_barcode '{rb}'")

    if errors:
        raise ValueError("Sample sheet errors:\n" + "\n".join(f"  {e}" for e in errors))

    return sample_ids


def resolve_seq(name, fw_dict, rv_dict, side):
    n = name.strip().upper()
    if side == "forward":
        return fw_dict.get(n, n)
    return rv_dict.get(n, n)


def main():
    args = parse_args()

    # Load protocol constants
    proto    = load_protocol(args.protocol_constants)
    fw_dict  = {k.upper(): v.upper() for k, v in proto["forward_barcodes"].items()
                if not k.startswith("_")}
    rv_dict  = {k.upper(): v.upper() for k, v in proto["reverse_barcodes"].items()
                if not k.startswith("_")}

    # Load sample sheet
    with open(args.samples, newline="") as f:
        rows = list(csv.DictReader(f))

    sample_ids = validate_sample_sheet(rows, fw_dict, rv_dict)
    print(f"[make_project] Loaded {len(rows)} samples from {args.samples}")

    # Create project directory
    proj_dir  = Path(args.project_dir).expanduser()
    proj_name = args.project_name or proj_dir.name

    # Pipeline subdirectory
    pipeline_dir = proj_dir / "pipeline"
    config_dir   = pipeline_dir / "config"
    scripts_dir  = pipeline_dir / "scripts"
    docs_dir     = pipeline_dir / "docs"
    examples_dir = pipeline_dir / "examples"
    tests_dir    = pipeline_dir / "tests"
    outputs_ex   = pipeline_dir / "outputs_example"
    fastq_dir    = proj_dir / "00_fastq"

    for d in [config_dir, scripts_dir, docs_dir, examples_dir, tests_dir,
              outputs_ex, fastq_dir]:
        d.mkdir(parents=True, exist_ok=True)
    print(f"[make_project] Created project structure under {proj_dir}")

    # Copy pipeline scripts and supporting files
    src_pipeline = Path(args.copy_pipeline)
    for sub in ["scripts", "tests", "docs", "examples"]:
        src_sub = src_pipeline / sub
        dst_sub = pipeline_dir / sub
        if src_sub.exists():
            for f in src_sub.iterdir():
                if f.is_file():
                    shutil.copy2(f, dst_sub / f.name)

    for f in ["README.md", "requirements.txt", "run_all.sh"]:
        src = src_pipeline / f
        if src.exists():
            shutil.copy2(src, pipeline_dir / f)

    # Copy static config files (except sample-specific ones)
    for f in ["hitom_protocol_constants.json"]:
        src = src_pipeline / "config" / f
        if src.exists():
            shutil.copy2(src, config_dir / f)
    print(f"[make_project] Copied pipeline files")

    # Generate sample_metadata.csv
    meta_cols   = ["sample_id", "editor", "target", "dpi", "replicate"]
    if "leaf" in rows[0]:
        meta_cols.append("leaf")
    meta_cols  += ["forward_barcode", "reverse_barcode"]

    with open(config_dir / "sample_metadata.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=meta_cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"[make_project] Wrote config/sample_metadata.csv ({len(rows)} samples)")

    # Generate barcode_map.csv (resolved sequences)
    bc_rows = []
    for row in rows:
        fw_seq = resolve_seq(row["forward_barcode"], fw_dict, rv_dict, "forward")
        rv_seq = resolve_seq(row["reverse_barcode"], fw_dict, rv_dict, "reverse")
        bc_rows.append({
            "sample_id":             row["sample_id"],
            "forward_barcode_name":  row["forward_barcode"],
            "forward_barcode_seq":   fw_seq,
            "reverse_barcode_name":  row["reverse_barcode"],
            "reverse_barcode_seq":   rv_seq,
        })
    with open(config_dir / "barcode_map.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(bc_rows[0].keys()))
        w.writeheader()
        w.writerows(bc_rows)
    print(f"[make_project] Wrote config/barcode_map.csv")

    # Generate targets.json template (one entry per unique target)
    unique_targets = list(dict.fromkeys(r["target"] for r in rows))
    targets_json   = dict(TARGETS_JSON_TEMPLATE)
    targets_json["targets"] = []
    for tgt in unique_targets:
        entry = dict(TARGET_ENTRY_TEMPLATE)
        entry["name"]           = f"{tgt}_MAIN"
        entry["applicable_to"]  = [tgt]
        entry["analysis_notes"] = f"Template entry for target {tgt}. Fill in motif sequences."
        targets_json["targets"].append(entry)

    with open(config_dir / "targets.json", "w") as f:
        json.dump(targets_json, f, indent=2)
    print(f"[make_project] Wrote config/targets.json template "
          f"({len(unique_targets)} target entries — FILL IN MOTIF SEQUENCES)")

    # Generate example run script
    fastq_path = args.raw_fastq_dir or str(fastq_dir)
    example_cmd = textwrap.dedent(f"""\
        #!/usr/bin/env bash
        # Example run for project: {proj_name}
        # Generated by make_new_hitom_project.py
        #
        # BEFORE RUNNING:
        #   1. Copy raw FASTQ.gz files to:  {fastq_path}
        #   2. Fill in motif sequences in:  {config_dir}/targets.json
        #   3. Install matplotlib/pandas:   pip install matplotlib pandas numpy
        #   4. Run this script

        set -euo pipefail
        PIPELINE_DIR="{pipeline_dir}"
        INPUT_DIR="{fastq_path}"
        OUTPUT_DIR="{proj_dir}/output"

        bash "$PIPELINE_DIR/run_all.sh" \\
            --input  "$INPUT_DIR" \\
            --output "$OUTPUT_DIR"

        echo "Done. Results in $OUTPUT_DIR"
    """)
    run_script = examples_dir / "run_this_project.sh"
    run_script.write_text(example_cmd)
    run_script.chmod(0o755)
    print(f"[make_project] Wrote examples/run_this_project.sh")

    # Summary report
    report = textwrap.dedent(f"""\
        # Project: {proj_name}
        # Generated: {__import__('datetime').date.today()}

        ## Quick start

        1. Copy raw FASTQs to:
               {fastq_path}

        2. Fill in motif sequences in:
               {config_dir}/targets.json
           (Search for FILL_IN_WT_MOTIF_SEQUENCE and replace with actual sequences)

        3. Run:
               bash {run_script}

        ## Sample summary
        Samples: {len(rows)}
        Targets: {', '.join(unique_targets)}
        Editors: {', '.join(sorted(set(r['editor'] for r in rows)))}
        DPI:     {', '.join(sorted(set(str(r['dpi']) for r in rows)))}

        ## Barcode usage
    """)
    for row in rows:
        report += f"  {row['sample_id']:<40s}  {row['forward_barcode']:>4s} × {row['reverse_barcode']:<4s}\n"

    report += textwrap.dedent(f"""
        ## Files generated
        {config_dir}/sample_metadata.csv   — sample sheet
        {config_dir}/barcode_map.csv       — resolved barcodes
        {config_dir}/targets.json          — motif definitions (FILL IN)
        {config_dir}/hitom_protocol_constants.json — HiTOM constants (do not modify)
        {examples_dir}/run_this_project.sh — example run command
    """)

    (proj_dir / "PROJECT_README.md").write_text(report)
    print(f"[make_project] Wrote PROJECT_README.md")

    print(f"\n[make_project] ✓ Project ready at {proj_dir}")
    print(f"[make_project]")
    print(f"[make_project] NEXT STEP: fill in motif sequences in")
    print(f"[make_project]   {config_dir}/targets.json")
    print(f"[make_project]")
    print(f"[make_project] Then run:")
    print(f"[make_project]   bash {run_script}")


if __name__ == "__main__":
    main()
