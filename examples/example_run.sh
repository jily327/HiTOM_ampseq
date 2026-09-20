#!/usr/bin/env bash
# =============================================================================
# example_run.sh  --  Minimal example to reproduce the ALSW analysis.
#
# Prerequisites:
#   - Python 3.8+
#   - Raw FASTQs copied to $HOME/projects/alsw/00_fastq/
#     (see README.md section "Input files")
#   - For figures: pip install matplotlib pandas
#
# Run from the project root:
#   bash pipeline/examples/example_run.sh
# =============================================================================

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

INPUT="$HOME/projects/alsw/00_fastq"
OUTPUT="$HOME/projects/alsw/rerun_output"

echo "Running ALSW pipeline..."
echo "  Input  : $INPUT"
echo "  Output : $OUTPUT"
echo ""

bash "$PIPELINE_DIR/run_all.sh" \
    --input  "$INPUT" \
    --output "$OUTPUT"

echo ""
echo "Done. Key result:"
echo ""
cat "$OUTPUT/04_edit/motif_count_summary.tsv" | \
    awk -F'\t' 'NR==1 || $1=="ALSW"' | \
    column -t
