#!/usr/bin/env python3
"""
09_run_crispresso2_batch.py  --  Run CRISPResso2 batch if installed.

*** OPTIONAL MODULE — NOT PART OF THE VALIDATED MOTIF-COUNT PIPELINE ***

Checks whether CRISPRessoBatch is available and runs it on the batch CSV
prepared by 08_prepare_crispresso2_inputs.py. If CRISPResso2 is not installed,
prints installation instructions and exits gracefully.

The main pipeline (run_all.sh) does NOT call this script by default.
Pass --run-crispresso2 to run_all.sh to enable.

Outputs (under <output>/09_crispresso2_results/):
  CRISPRessoBatch output directories (one per sample/target)
  crispresso2_run.log
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils

CRISPRESSO_CMDS = ["CRISPRessoBatch", "crispressoBatch", "crispresso2_batch"]


def find_crispresso():
    """Return the CRISPRessoBatch executable path, or None."""
    for cmd in CRISPRESSO_CMDS:
        found = shutil.which(cmd)
        if found:
            return found
    return None


def print_install_instructions():
    print("""
[09_crispresso2] CRISPRessoBatch not found in PATH.

To install CRISPResso2:

  Option A — conda (recommended):
    conda create -n crispresso2 -c conda-forge -c bioconda crispresso2 -y
    conda activate crispresso2

  Option B — pip:
    pip install CRISPResso2

  Option C — Docker:
    docker pull pinellolab/crispresso2
    docker run pinellolab/crispresso2 CRISPRessoBatch --help

After installation, re-run this script or run manually:
    CRISPRessoBatch --batch_settings <output>/08_crispresso2_inputs/CRISPRessoBatch_input.csv

See pipeline/docs/crispresso2_local_usage.md for full guidance.
""")


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input",    required=True,
                   help="Root output dir (expects <input>/08_crispresso2_inputs/)")
    p.add_argument("--output",   required=True,
                   help="Root output dir (results go to <output>/09_crispresso2_results/)")
    p.add_argument("--threads",  type=int, default=4,
                   help="CRISPRessoBatch --n_processes [default: 4]")
    p.add_argument("--min-freq", type=float, default=0.05, dest="min_freq",
                   help="--min_frequency_alleles_around_cut_to_plot [default: 0.05]")
    p.add_argument("--extra-args", default="", dest="extra_args",
                   help="Additional CRISPRessoBatch arguments (quoted string)")
    return p.parse_args()


def main():
    args       = parse_args()
    batch_csv  = Path(args.input) / "08_crispresso2_inputs" / "CRISPRessoBatch_input.csv"
    out_dir    = Path(args.output) / "09_crispresso2_results"
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path   = out_dir / "crispresso2_run.log"

    # Check CRISPResso2 is installed FIRST — exit 0 gracefully if absent
    crispresso_exe = find_crispresso()
    if crispresso_exe is None:
        print_install_instructions()
        # Graceful exit — do not cause run_all.sh to fail
        sys.exit(0)

    # CRISPResso2 is present — now validate inputs
    if not batch_csv.exists():
        print(f"[09_crispresso2] ERROR: batch CSV not found: {batch_csv}",
              file=sys.stderr)
        print(f"[09_crispresso2]   Run 08_prepare_crispresso2_inputs.py first.",
              file=sys.stderr)
        sys.exit(1)

    # Check for placeholder amplicon sequences
    rows = utils.read_csv_file(batch_csv)
    n_placeholder = sum(1 for r in rows
                        if "FILL_IN_AMPLICON_SEQUENCE" in r.get("amplicon_seq", ""))
    if n_placeholder > 0:
        print(f"[09_crispresso2] ERROR: {n_placeholder} rows have "
              f"FILL_IN_AMPLICON_SEQUENCE in amplicon_seq.", file=sys.stderr)
        print(f"[09_crispresso2]   Edit {batch_csv} and fill in amplicon "
              f"sequences before running.", file=sys.stderr)
        sys.exit(1)

    print(f"[09_crispresso2] Found CRISPRessoBatch: {crispresso_exe}")
    print(f"[09_crispresso2] Input CSV : {batch_csv}")
    print(f"[09_crispresso2] Output dir: {out_dir}")
    print(f"[09_crispresso2] Threads   : {args.threads}")

    cmd = [
        crispresso_exe,
        "--batch_settings", str(batch_csv),
        "--output_folder",  str(out_dir),
        "--n_processes",    str(args.threads),
        "--min_frequency_alleles_around_cut_to_plot", str(args.min_freq),
    ]
    if args.extra_args:
        cmd.extend(args.extra_args.split())

    print(f"[09_crispresso2] Command: {' '.join(cmd)}")

    with open(log_path, "w", encoding=utils.TEXT_ENCODING) as log:
        log.write(f"Command: {' '.join(cmd)}\n\n")
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)

    if result.returncode == 0:
        print(f"[09_crispresso2] CRISPRessoBatch completed successfully.")
        print(f"[09_crispresso2] Results in {out_dir}")
        print(f"[09_crispresso2] Log: {log_path}")
    else:
        print(f"[09_crispresso2] ERROR: CRISPRessoBatch exited with code "
              f"{result.returncode}", file=sys.stderr)
        print(f"[09_crispresso2] See log: {log_path}", file=sys.stderr)
        sys.exit(result.returncode)


if __name__ == "__main__":
    main()
