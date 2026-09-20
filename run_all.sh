#!/usr/bin/env bash
# =============================================================================
# run_all.sh  --  ALSW HiTOM amplicon-seq pipeline orchestrator
#
# Usage:
#   bash pipeline/run_all.sh \
#       --input  /path/to/raw_fastq_dir \
#       --output /path/to/output_dir
#
# Optional config overrides:
#   --metadata         /custom/sample_metadata.csv
#   --config           /custom/targets.json
#   --protocol-constants /custom/hitom_protocol_constants.json
#
# All default config files live in pipeline/config/.
# No hardcoded absolute paths are required for normal use.
#
# Steps (validated, always run):
#   01  Barcode demux + QC
#   02  Primer-dimer / short-read QC (informational)
#   03  Allele assignment
#   04  Edit classification → motif_count_summary.tsv
#   05  Summary tables
# Optional steps (run only if flag is passed):
#   07  Alignment-based edit classification     (--run-alignment)
#   08  Prepare CRISPResso2 inputs              (--prepare-crispresso2)
#   09  Run CRISPResso2 batch                   (--run-crispresso2)
#   10  Motif vs alignment comparison           (--run-alignment)
#   11  Allele-aware precise edit               (--run-allele-aware)
# Always runs last (benefits from any optional outputs produced above):
#   06  Figures
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$SCRIPT_DIR/scripts"

# ---------------------------------------------------------------------------
# Defaults — point to pipeline/config/
# ---------------------------------------------------------------------------
DEFAULT_METADATA="$SCRIPT_DIR/config/sample_metadata.csv"
DEFAULT_CONFIG="$SCRIPT_DIR/config/targets.json"
DEFAULT_PROTOCOL="$SCRIPT_DIR/config/hitom_protocol_constants.json"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
INPUT=""
OUTPUT=""
METADATA="$DEFAULT_METADATA"
CONFIG="$DEFAULT_CONFIG"
PROTOCOL="$DEFAULT_PROTOCOL"
# Optional module flags (off by default)
RUN_ALIGNMENT=false
PREPARE_CRISPRESSO2=false
RUN_CRISPRESSO2=false
RUN_ALLELE_AWARE=false

usage() {
    sed -n '3,15p' "$0"
    echo ""
    echo "Optional module flags (not run by default):"
    echo "  --run-alignment        Run 07_alignment_edit_classifier.py"
    echo "  --prepare-crispresso2  Run 08_prepare_crispresso2_inputs.py"
    echo "  --run-crispresso2      Run 09_run_crispresso2_batch.py (requires CRISPResso2)"
    echo "  --run-allele-aware     Run 11_allele_aware_precise_edit.py"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)               INPUT="$2";    shift 2 ;;
        --output)              OUTPUT="$2";   shift 2 ;;
        --metadata)            METADATA="$2"; shift 2 ;;
        --config)              CONFIG="$2";   shift 2 ;;
        --protocol-constants)  PROTOCOL="$2"; shift 2 ;;
        --run-alignment)       RUN_ALIGNMENT=true; shift ;;
        --prepare-crispresso2) PREPARE_CRISPRESSO2=true; shift ;;
        --run-crispresso2)     RUN_CRISPRESSO2=true; PREPARE_CRISPRESSO2=true; shift ;;
        --run-allele-aware)    RUN_ALLELE_AWARE=true; shift ;;
        -h|--help)             usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo "ERROR: --input and --output are required."
    usage
fi

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
if [[ ! -d "$INPUT" ]]; then
    echo "ERROR: Input directory not found: $INPUT"
    exit 1
fi
if [[ ! -f "$METADATA" ]]; then
    echo "ERROR: Metadata file not found: $METADATA"
    exit 1
fi
if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi
if [[ ! -f "$PROTOCOL" ]]; then
    echo "ERROR: Protocol constants file not found: $PROTOCOL"
    exit 1
fi

mkdir -p "$OUTPUT"

echo "======================================================================"
echo "ALSW HiTOM Amplicon-seq Pipeline"
echo "======================================================================"
echo "  Input FASTQ dir  : $INPUT"
echo "  Output dir       : $OUTPUT"
echo "  Metadata         : $METADATA"
echo "  Targets config   : $CONFIG"
echo "  Protocol consts  : $PROTOCOL"
echo "  Started          : $(date '+%Y-%m-%d %H:%M:%S')"
echo "======================================================================"

PYTHON=$(command -v python3 || command -v python)
if [[ -z "$PYTHON" ]]; then
    echo "ERROR: python3 not found in PATH."
    exit 1
fi
echo "  Python           : $($PYTHON --version 2>&1)"
echo ""

# ---------------------------------------------------------------------------
# Step 01  Barcode demux + QC
# ---------------------------------------------------------------------------
echo "--- Step 01: Barcode demux + QC ---"
"$PYTHON" "$SCRIPTS_DIR/01_demux_and_qc.py" \
    --input               "$INPUT" \
    --output              "$OUTPUT" \
    --metadata            "$METADATA" \
    --protocol-constants  "$PROTOCOL"
echo ""

# ---------------------------------------------------------------------------
# Step 02  Primer-dimer / short-read QC  (informational only)
# ---------------------------------------------------------------------------
echo "--- Step 02: Primer-dimer QC ---"
"$PYTHON" "$SCRIPTS_DIR/02_primer_dimer_filter.py" \
    --input  "$OUTPUT" \
    --output "$OUTPUT" \
    --metadata "$METADATA"
echo ""

# ---------------------------------------------------------------------------
# Step 03  Allele assignment
# ---------------------------------------------------------------------------
echo "--- Step 03: Allele assignment ---"
"$PYTHON" "$SCRIPTS_DIR/03_allele_assignment.py" \
    --input  "$OUTPUT" \
    --output "$OUTPUT" \
    --metadata "$METADATA" \
    --config   "$CONFIG"
echo ""

# ---------------------------------------------------------------------------
# Step 04  Edit classification  →  motif_count_summary.tsv
# ---------------------------------------------------------------------------
echo "--- Step 04: Edit classification ---"
"$PYTHON" "$SCRIPTS_DIR/04_edit_classification.py" \
    --input  "$OUTPUT" \
    --output "$OUTPUT" \
    --metadata "$METADATA" \
    --config   "$CONFIG"
echo ""

# ---------------------------------------------------------------------------
# Step 05  Summary tables
# ---------------------------------------------------------------------------
echo "--- Step 05: Summary tables ---"
"$PYTHON" "$SCRIPTS_DIR/05_make_summary_tables.py" \
    --input  "$OUTPUT" \
    --output "$OUTPUT" \
    --metadata "$METADATA"
echo ""

# ---------------------------------------------------------------------------
# Optional Step 07  Alignment-based edit classification
# ---------------------------------------------------------------------------
if [[ "$RUN_ALIGNMENT" == "true" ]]; then
    echo "--- Step 07: Alignment-based edit classification (OPTIONAL) ---"
    echo "    *** NOT VALIDATED — exploratory analysis only ***"
    "$PYTHON" "$SCRIPTS_DIR/07_alignment_edit_classifier.py" \
        --input    "$OUTPUT" \
        --output   "$OUTPUT" \
        --metadata "$METADATA" \
        --config   "$CONFIG"
    echo ""
fi

# ---------------------------------------------------------------------------
# Optional Step 08  Prepare CRISPResso2 inputs
# ---------------------------------------------------------------------------
if [[ "$PREPARE_CRISPRESSO2" == "true" ]]; then
    echo "--- Step 08: Prepare CRISPResso2 inputs (OPTIONAL) ---"
    "$PYTHON" "$SCRIPTS_DIR/08_prepare_crispresso2_inputs.py" \
        --input    "$OUTPUT" \
        --output   "$OUTPUT" \
        --metadata "$METADATA" \
        --config   "$CONFIG"
    echo ""
fi

# ---------------------------------------------------------------------------
# Optional Step 09  Run CRISPResso2 (requires installation)
# ---------------------------------------------------------------------------
if [[ "$RUN_CRISPRESSO2" == "true" ]]; then
    echo "--- Step 09: Run CRISPResso2 batch (OPTIONAL) ---"
    "$PYTHON" "$SCRIPTS_DIR/09_run_crispresso2_batch.py" \
        --input  "$OUTPUT" \
        --output "$OUTPUT"
    echo ""
fi

# ---------------------------------------------------------------------------
# Optional Step 10  Motif vs alignment comparison (runs if both 04 + 07 done)
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Optional Step 11  Allele-aware precise edit quantification
# ---------------------------------------------------------------------------
if [[ "$RUN_ALLELE_AWARE" == "true" ]]; then
    echo "--- Step 11: Allele-aware precise edit quantification (OPTIONAL) ---"
    echo "    *** NOT VALIDATED — exploratory analysis only ***"
    "$PYTHON" "$SCRIPTS_DIR/11_allele_aware_precise_edit.py" \
        --input    "$OUTPUT" \
        --output   "$OUTPUT" \
        --metadata "$METADATA" \
        --config   "$CONFIG"
    echo ""
fi

if [[ "$RUN_ALIGNMENT" == "true" ]]; then
    echo "--- Step 10: Motif vs alignment comparison (OPTIONAL) ---"
    "$PYTHON" "$SCRIPTS_DIR/10_compare_motif_vs_alignment.py" \
        --input    "$OUTPUT" \
        --output   "$OUTPUT" \
        --metadata "$METADATA"
    echo ""
fi

# ---------------------------------------------------------------------------
# Step 06  Figures — runs LAST so all optional outputs (incl. allele-aware)
#           are already present and fig8 is generated when available.
#           (requires matplotlib + pandas + numpy; skipped with warning if absent)
# ---------------------------------------------------------------------------
echo "--- Step 06: Figures ---"
if "$PYTHON" -c "import matplotlib, pandas, numpy" 2>/dev/null; then
    "$PYTHON" "$SCRIPTS_DIR/06_make_figures.py" \
        --input    "$OUTPUT" \
        --output   "$OUTPUT" \
        --metadata "$METADATA" \
        --format   "pdf,png,svg" \
        --dpi      300
else
    echo "[06_figures] SKIPPED: matplotlib, pandas, and/or numpy not installed."
    echo "  Install with:  pip install matplotlib pandas numpy"
    echo "  Then run step 06 manually:"
    echo "    python3 $SCRIPTS_DIR/06_make_figures.py \\"
    echo "        --input $OUTPUT --output $OUTPUT --format pdf,png,svg --dpi 300"
fi
echo ""

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "======================================================================"
echo "Pipeline complete"
echo "  Finished         : $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo "Key outputs:"
echo "  motif_count_summary.tsv : $OUTPUT/04_edit/motif_count_summary.tsv"
echo "  final_sample_summary    : $OUTPUT/05_summary/final_sample_summary.csv"
echo "  Figures                 : $OUTPUT/06_figures/"
if [[ "$RUN_ALIGNMENT" == "true" ]]; then
echo "  Alignment summary       : $OUTPUT/07_alignment/alignment_sample_summary.csv"
echo "  Motif vs alignment      : $OUTPUT/10_comparison/motif_vs_alignment_comparison.csv"
fi
if [[ "$RUN_ALLELE_AWARE" == "true" ]]; then
echo "  Allele-aware summary    : $OUTPUT/11_allele_aware/allele_aware_precise_edit_summary.csv"
fi
echo "======================================================================"
