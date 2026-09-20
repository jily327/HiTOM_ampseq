#!/usr/bin/env bash
# =============================================================================
# smoke_test.sh  --  Pipeline smoke test and reproducibility check.
#
# Runs the full pipeline into a fresh output directory and compares
# regenerated outputs against gold-standard reference files in
# pipeline/outputs_example/.
#
# Usage:
#   bash pipeline/tests/smoke_test.sh
#
# Exit codes:
#   0 = all checks passed
#   1 = one or more checks failed
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GOLD_DIR="$SCRIPT_DIR/outputs_example"
INPUT_DIR="$HOME/projects/alsw/00_fastq"
SMOKE_OUT="$HOME/projects/alsw/smoke_test_output_$(date +%Y%m%d_%H%M%S)"

RED='\033[0;31m'
GRN='\033[0;32m'
YLW='\033[0;33m'
NC='\033[0m'

PASS=0
FAIL=0
WARN=0

pass() { echo -e "${GRN}PASS${NC}  $1"; ((PASS++)) || true; }
fail() { echo -e "${RED}FAIL${NC}  $1"; ((FAIL++)) || true; }
warn() { echo -e "${YLW}WARN${NC}  $1"; ((WARN++)) || true; }

echo "======================================================================"
echo "ALSW Portable Pipeline  —  Smoke Test"
echo "  Pipeline dir : $SCRIPT_DIR"
echo "  Input FASTQs : $INPUT_DIR"
echo "  Output dir   : $SMOKE_OUT"
echo "  Gold standard: $GOLD_DIR"
echo "  Started      : $(date '+%Y-%m-%d %H:%M:%S')"
echo "======================================================================"
echo ""

# ---------------------------------------------------------------------------
# 0. Preflight checks
# ---------------------------------------------------------------------------
echo "--- Preflight ---"

if [[ ! -d "$INPUT_DIR" ]]; then
    fail "Input directory not found: $INPUT_DIR"
    echo "  Ensure raw FASTQs are in $INPUT_DIR"
    exit 1
fi
pass "Input directory exists: $INPUT_DIR"

R1_COUNT=$(ls "$INPUT_DIR"/*_R1*.fastq.gz 2>/dev/null | wc -l)
R2_COUNT=$(ls "$INPUT_DIR"/*_R2*.fastq.gz 2>/dev/null | wc -l)
if [[ "$R1_COUNT" -eq 1 && "$R2_COUNT" -eq 1 ]]; then
    pass "Exactly 1 R1 and 1 R2 FASTQ.gz found"
else
    fail "Expected 1 R1 + 1 R2 FASTQ.gz, found R1=$R1_COUNT R2=$R2_COUNT"
fi

if [[ ! -f "$GOLD_DIR/motif_count_summary.tsv" ]]; then
    fail "Gold-standard motif_count_summary.tsv not found in $GOLD_DIR"
    exit 1
fi
pass "Gold-standard motif_count_summary.tsv found"

GOLD_ROWS=$(wc -l < "$GOLD_DIR/motif_count_summary.tsv")
pass "Gold standard has $GOLD_ROWS lines (including header)"
echo ""

# ---------------------------------------------------------------------------
# 1. Check for hardcoded absolute paths in scripts
# ---------------------------------------------------------------------------
echo "--- Hardcoded absolute path check ---"
SCRIPTS_DIR="$SCRIPT_DIR/scripts"
HARDCODED=0
while IFS= read -r -d '' pyfile; do
    # Grep for /home/ /mnt/ /Users/ /root/ patterns that are NOT in comments
    HITS=$(grep -n '["'"'"']\(/home\|/mnt\|/Users\|/root\|/tmp\)' "$pyfile" 2>/dev/null \
           | grep -v '^\s*#' || true)
    if [[ -n "$HITS" ]]; then
        fail "Hardcoded path in $(basename "$pyfile"):"
        echo "$HITS" | sed 's/^/    /'
        ((HARDCODED++)) || true
    fi
done < <(find "$SCRIPTS_DIR" -name "*.py" -print0)

# Also check run_all.sh
HITS=$(grep -n '["'"'"']\(/home\|/mnt\|/Users\|/root\)' "$SCRIPT_DIR/run_all.sh" 2>/dev/null \
       | grep -v '^\s*#' || true)
if [[ -n "$HITS" ]]; then
    fail "Hardcoded path in run_all.sh:"
    echo "$HITS" | sed 's/^/    /'
    ((HARDCODED++)) || true
fi

if [[ "$HARDCODED" -eq 0 ]]; then
    pass "No hardcoded absolute paths found in scripts"
fi
echo ""

# ---------------------------------------------------------------------------
# 2. Run the full pipeline
# ---------------------------------------------------------------------------
echo "--- Running pipeline ---"
mkdir -p "$SMOKE_OUT"

if bash "$SCRIPT_DIR/run_all.sh" \
        --input  "$INPUT_DIR" \
        --output "$SMOKE_OUT" 2>&1 | tee "$SMOKE_OUT/pipeline_run.log"; then
    pass "Pipeline exited cleanly"
else
    fail "Pipeline returned non-zero exit code"
    echo "  See log: $SMOKE_OUT/pipeline_run.log"
    exit 1
fi
echo ""

# ---------------------------------------------------------------------------
# 3. Check expected output files exist
# ---------------------------------------------------------------------------
echo "--- Output file existence ---"

EXPECTED_FILES=(
    "01_demux/demux_counts.tsv"
    "01_demux/tail_trim_summary.tsv"
    "01_demux/barcode_assignment_summary.csv"
    "01_demux/read_length_summary.csv"
    "01_demux/barcode_map_resolved.csv"
    "01_demux/md5_summary.txt"
    "02_qc/primer_dimer_qc_summary.csv"
    "02_qc/per_sample_length_histogram.csv"
    "03_allele/allele_counts.tsv"
    "04_edit/motif_count_summary.tsv"
    "04_edit/edit_counts.tsv"
    "05_summary/final_sample_summary.csv"
    "05_summary/per_condition_summary.csv"
    "05_summary/wt_background_summary.csv"
    "05_summary/allele_specific_summary.csv"
    "05_summary/sensitivity_denominators.csv"
    "05_summary/primer_dimer_qc_summary.csv"
    "06_figures/figure_legends.md"
)

for f in "${EXPECTED_FILES[@]}"; do
    if [[ -f "$SMOKE_OUT/$f" ]]; then
        pass "Exists: $f"
    else
        fail "Missing: $f"
    fi
done
echo ""

# ---------------------------------------------------------------------------
# 4. Compare motif_count_summary.tsv against gold standard
# ---------------------------------------------------------------------------
echo "--- motif_count_summary.tsv numerical comparison ---"
REGEN="$SMOKE_OUT/04_edit/motif_count_summary.tsv"
GOLD="$GOLD_DIR/motif_count_summary.tsv"

if [[ ! -f "$REGEN" ]]; then
    fail "Regenerated motif_count_summary.tsv not found; cannot compare"
else
    # Row count
    REGEN_ROWS=$(wc -l < "$REGEN")
    if [[ "$REGEN_ROWS" -eq "$GOLD_ROWS" ]]; then
        pass "Row count matches: $REGEN_ROWS rows"
    else
        fail "Row count mismatch: regenerated=$REGEN_ROWS gold=$GOLD_ROWS"
    fi

    # Exact diff
    if diff -q "$GOLD" "$REGEN" > /dev/null 2>&1; then
        pass "motif_count_summary.tsv is byte-identical to gold standard"
    else
        # Show differing lines
        DIFF_LINES=$(diff "$GOLD" "$REGEN" | grep "^[<>]" | wc -l)
        fail "motif_count_summary.tsv differs from gold standard ($DIFF_LINES differing lines)"
        echo ""
        echo "  First differences (gold < vs regen >):"
        diff "$GOLD" "$REGEN" | grep "^[<>]" | head -20 | sed 's/^/    /'

        # Compare numeric columns for ALSW (key biological result)
        echo ""
        echo "  ALSW rows in gold:"
        grep "^ALSW" "$GOLD" | awk -F'\t' '{printf "  %-35s  desired=%s  pct=%s\n", $2, $5, $9}'
        echo "  ALSW rows in regen:"
        grep "^ALSW" "$REGEN" | awk -F'\t' '{printf "  %-35s  desired=%s  pct=%s\n", $2, $5, $9}'
    fi
fi
echo ""

# ---------------------------------------------------------------------------
# 5. Spot-check key ALSW values against known-good results
# ---------------------------------------------------------------------------
echo "--- ALSW key value spot-check ---"
# Known values from analysis_summary_updated.txt:
#   ePPEmax_ALSW_6dpi_pool: ~23.6339%
#   PE6c_ALSW_6dpi_pool:    ~3.5203%
#   WT_ALSW_6dpi_pool:       ~0.3712% (background)

check_value() {
    local label="$1" file="$2" grep_pat="$3" col="$4" expected="$5"
    local val
    val=$(grep "$grep_pat" "$file" 2>/dev/null | awk -F'\t' "{print \$$col}" | head -1)
    if [[ "$val" == "$expected" ]]; then
        pass "$label = $val"
    else
        fail "$label expected=$expected got=$val"
    fi
}

if [[ -f "$REGEN" ]]; then
    check_value "ePPEmax_ALSW_6dpi desired_pct" \
        "$REGEN" "ePPEmax_ALSW_6dpi_pool" 9 "23.6339"
    check_value "PE6c_ALSW_6dpi desired_pct" \
        "$REGEN" "PE6c_ALSW_6dpi_pool" 9 "3.5203"
    check_value "WT_ALSW_6dpi desired_pct (background)" \
        "$REGEN" "WT_ALSW_6dpi_pool" 9 "0.3712"
    check_value "ePPEmax_ALSW_3dpi desired_pct" \
        "$REGEN" "ePPEmax_ALSW_3dpi_pool" 9 "13.9483"
fi
echo ""

# ---------------------------------------------------------------------------
# 6. Figure check (PDF, SVG, PNG subdirectories + source_data)
# ---------------------------------------------------------------------------
echo "--- Figure check ---"
FIG_DIR="$SMOKE_OUT/06_figures"
PYTHON=$(command -v python3 || command -v python || echo "")
if [[ -n "$PYTHON" ]] && "$PYTHON" -c "import matplotlib, pandas, numpy" 2>/dev/null; then
    PDF_COUNT=$(ls "$FIG_DIR"/pdf/*.pdf 2>/dev/null | wc -l)
    PNG_COUNT=$(ls "$FIG_DIR"/png/*.png 2>/dev/null | wc -l)
    SVG_COUNT=$(ls "$FIG_DIR"/svg/*.svg 2>/dev/null | wc -l)
    SRC_COUNT=$(ls "$FIG_DIR"/source_data/*.csv 2>/dev/null | wc -l)

    if [[ "$PDF_COUNT" -ge 7 ]]; then
        pass "$PDF_COUNT PDF figures in 06_figures/pdf/"
    else
        fail "Expected ≥7 PDFs, found $PDF_COUNT in $FIG_DIR/pdf/"
    fi

    if [[ "$PNG_COUNT" -ge 7 ]]; then
        pass "$PNG_COUNT PNG figures in 06_figures/png/"
    else
        fail "Expected ≥7 PNGs, found $PNG_COUNT in $FIG_DIR/png/"
    fi

    if [[ "$SVG_COUNT" -ge 7 ]]; then
        pass "$SVG_COUNT SVG figures in 06_figures/svg/"
    else
        fail "Expected ≥7 SVGs, found $SVG_COUNT in $FIG_DIR/svg/"
    fi

    if [[ "$SRC_COUNT" -ge 6 ]]; then
        pass "$SRC_COUNT source-data CSVs in 06_figures/source_data/"
    else
        fail "Expected ≥6 source CSVs, found $SRC_COUNT in $FIG_DIR/source_data/"
    fi

    if [[ -f "$FIG_DIR/figure_legends.md" ]]; then
        pass "figure_legends.md exists"
    else
        fail "figure_legends.md not found in $FIG_DIR"
    fi

    echo "  Figure stems:"
    for f in "$FIG_DIR"/pdf/*.pdf; do
        [[ -f "$f" ]] && echo "    $(basename "${f%.pdf}")"
    done
else
    warn "matplotlib/pandas/numpy not installed — figures not tested"
    warn "  Run: pip install matplotlib pandas numpy"
fi
echo ""

# ---------------------------------------------------------------------------
# 7. Summary
# ---------------------------------------------------------------------------
echo "======================================================================"
echo "Smoke test summary"
echo "  Finished : $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Output   : $SMOKE_OUT"
echo ""
echo -e "  ${GRN}PASS${NC} : $PASS"
echo -e "  ${RED}FAIL${NC} : $FAIL"
echo -e "  ${YLW}WARN${NC} : $WARN"
echo ""

if [[ "$FAIL" -eq 0 ]]; then
    echo -e "${GRN}All checks passed.${NC}"
    echo ""
    echo "To reproduce the analysis at any time:"
    echo ""
    echo "  bash pipeline/run_all.sh \\"
    echo "      --input  \"\$HOME/projects/alsw/00_fastq\" \\"
    echo "      --output \"\$HOME/projects/alsw/rerun_output\""
    echo ""
    exit 0
else
    echo -e "${RED}$FAIL check(s) FAILED. See output above.${NC}"
    exit 1
fi
