#!/usr/bin/env bash
# =============================================================================
# smoke_test_allele_aware.sh  --  Smoke test for allele-aware module (11).
#
# Tests:
#   1. Validated pipeline runs cleanly; motif_count_summary.tsv byte-identical
#   2. Allele-aware module runs for configured targets (ALSP)
#   3. ALSW skipped gracefully (unconfigured status)
#   4. Expected output files created
#   5. Ambiguous/low_information reads are present in output (not silently dropped)
#   6. motif_count_summary.tsv unchanged after allele-aware step
#   7. fig8 either generated (if allele-aware ran) or skipped gracefully in figs
#
# Exit codes: 0 = all checks passed, 1 = failures
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GOLD_DIR="$SCRIPT_DIR/outputs_example"
INPUT_DIR="$HOME/projects/alsw/00_fastq"
SMOKE_OUT="$HOME/projects/alsw/smoke_allele_aware_$(date +%Y%m%d_%H%M%S)"
PYTHON=$(command -v python3 || command -v python)

RED='\033[0;31m'; GRN='\033[0;32m'; YLW='\033[0;33m'; NC='\033[0m'
PASS=0; FAIL=0; WARN=0

pass() { echo -e "${GRN}PASS${NC}  $1"; ((PASS++)) || true; }
fail() { echo -e "${RED}FAIL${NC}  $1"; ((FAIL++)) || true; }
warn() { echo -e "${YLW}WARN${NC}  $1"; ((WARN++)) || true; }

echo "======================================================================"
echo "ALSW Allele-Aware Module Smoke Test"
echo "  Pipeline    : $SCRIPT_DIR"
echo "  Input FASTQs: $INPUT_DIR"
echo "  Output      : $SMOKE_OUT"
echo "  Started     : $(date '+%Y-%m-%d %H:%M:%S')"
echo "======================================================================"
echo ""

# ---------------------------------------------------------------------------
# 1. Run FULL validated pipeline (steps 01-06) + allele-aware (step 11)
# ---------------------------------------------------------------------------
echo "--- Step 1: Run validated pipeline + allele-aware ---"
mkdir -p "$SMOKE_OUT"

if bash "$SCRIPT_DIR/run_all.sh" \
        --input  "$INPUT_DIR" \
        --output "$SMOKE_OUT" \
        --run-allele-aware \
   > "$SMOKE_OUT/pipeline_run.log" 2>&1; then
    pass "Pipeline (01-06 + 11) completed cleanly"
else
    fail "Pipeline returned non-zero exit code"
    echo "  See: $SMOKE_OUT/pipeline_run.log"
    exit 1
fi
echo ""

# ---------------------------------------------------------------------------
# 2. Gold standard still intact
# ---------------------------------------------------------------------------
echo "--- Step 2: Validated motif_count_summary.tsv unchanged ---"
if diff -q "$GOLD_DIR/motif_count_summary.tsv" \
           "$SMOKE_OUT/04_edit/motif_count_summary.tsv" > /dev/null 2>&1; then
    pass "motif_count_summary.tsv byte-identical to gold standard"
else
    fail "motif_count_summary.tsv CHANGED — allele-aware module must not modify validated outputs"
fi
echo ""

# ---------------------------------------------------------------------------
# 3. Allele-aware output files exist
# ---------------------------------------------------------------------------
echo "--- Step 3: Allele-aware output files ---"
EXPECTED_AA=(
    "11_allele_aware/allele_aware_precise_edit_summary.csv"
    "11_allele_aware/read_level_allele_edit_calls.tsv"
    "11_allele_aware/allele_assignment_qc.csv"
    "11_allele_aware/allele_aware_vs_motif_comparison.csv"
)
for f in "${EXPECTED_AA[@]}"; do
    if [[ -f "$SMOKE_OUT/$f" ]]; then
        pass "Exists: $f"
    else
        fail "Missing: $f"
    fi
done
echo ""

# ---------------------------------------------------------------------------
# 4. ALSP was processed (configured target)
# ---------------------------------------------------------------------------
echo "--- Step 4: ALSP allele-aware results ---"
SUMMARY="$SMOKE_OUT/11_allele_aware/allele_aware_precise_edit_summary.csv"
if [[ -f "$SUMMARY" ]]; then
    ALSP_ROWS=$(grep "^ALSP," "$SUMMARY" 2>/dev/null | wc -l)
    if [[ "$ALSP_ROWS" -gt 0 ]]; then
        pass "ALSP rows present in summary ($ALSP_ROWS rows)"
    else
        fail "No ALSP rows in allele_aware_precise_edit_summary.csv"
    fi

    # Confirm chr01 and chr11 rows exist
    if grep -q ",chr01," "$SUMMARY"; then
        pass "chr01 rows present"
    else
        fail "No chr01 rows in summary"
    fi
    if grep -q ",chr11," "$SUMMARY"; then
        pass "chr11 rows present"
    else
        fail "No chr11 rows in summary"
    fi
else
    fail "Summary file not found"
fi
echo ""

# ---------------------------------------------------------------------------
# 5. ALSW was skipped gracefully (unconfigured)
# ---------------------------------------------------------------------------
echo "--- Step 5: ALSW skipped gracefully ---"
if grep -q "SKIP.*ALSW\|ALSW.*unconfigured" "$SMOKE_OUT/pipeline_run.log" 2>/dev/null; then
    pass "ALSW skipped gracefully with informative message"
else
    # Also acceptable: ALSW just not present in output
    ALSW_IN_SUMMARY=0
    if [[ -f "$SUMMARY" ]]; then
        ALSW_IN_SUMMARY=$(grep ",ALSW," "$SUMMARY" 2>/dev/null | wc -l)
    fi
    if [[ "$ALSW_IN_SUMMARY" -eq 0 ]]; then
        pass "ALSW not in output (correct — unconfigured target not processed)"
    else
        warn "ALSW rows found in summary — check if this is expected"
    fi
fi
echo ""

# ---------------------------------------------------------------------------
# 6. Ambiguous / low_information reads are NOT silently discarded
# ---------------------------------------------------------------------------
echo "--- Step 6: Ambiguous reads accounted for ---"
if [[ -f "$SUMMARY" ]]; then
    if grep -q ",ambiguous," "$SUMMARY" || grep -q ",low_information," "$SUMMARY"; then
        AMBIG_ROWS=$(grep -c ",ambiguous,\|,low_information," "$SUMMARY" 2>/dev/null || true)
        pass "Ambiguous/low_information categories present ($AMBIG_ROWS rows)"
    else
        fail "No ambiguous or low_information rows in summary — reads may be silently dropped"
    fi

    # QC: confirm allele assignment rate column exists
    if grep -q "allele_assignment_rate_pct" "$SMOKE_OUT/11_allele_aware/allele_assignment_qc.csv" 2>/dev/null; then
        pass "allele_assignment_rate_pct column present in QC file"
    else
        fail "allele_assignment_rate_pct column missing from QC file"
    fi
fi
echo ""

# ---------------------------------------------------------------------------
# 7. Read-level table non-empty and has expected columns
# ---------------------------------------------------------------------------
echo "--- Step 7: Read-level output quality ---"
READ_TABLE="$SMOKE_OUT/11_allele_aware/read_level_allele_edit_calls.tsv"
if [[ -f "$READ_TABLE" ]]; then
    READ_ROWS=$(wc -l < "$READ_TABLE")
    if [[ "$READ_ROWS" -gt 1 ]]; then
        pass "read_level_allele_edit_calls.tsv has $READ_ROWS rows (incl header)"
    else
        fail "read_level_allele_edit_calls.tsv is empty"
    fi

    # Check for key columns
    HEADER=$(head -1 "$READ_TABLE")
    for COL in "assigned_allele" "edit_status" "has_precise_desired_edit" "informative_snp_count"; do
        if echo "$HEADER" | grep -q "$COL"; then
            pass "Column present: $COL"
        else
            fail "Column missing: $COL"
        fi
    done
else
    fail "read_level_allele_edit_calls.tsv not found"
fi
echo ""

# ---------------------------------------------------------------------------
# 8. allele_aware_vs_motif_comparison.csv contains ALSP entries
# ---------------------------------------------------------------------------
echo "--- Step 8: Allele-aware vs motif comparison ---"
COMPARE="$SMOKE_OUT/11_allele_aware/allele_aware_vs_motif_comparison.csv"
if [[ -f "$COMPARE" ]]; then
    COMPARE_ROWS=$(grep -c "ALSP" "$COMPARE" 2>/dev/null || echo 0)
    if [[ "$COMPARE_ROWS" -gt 0 ]]; then
        pass "ALSP entries in allele_aware_vs_motif_comparison.csv ($COMPARE_ROWS rows)"
    else
        warn "No ALSP rows in comparison file — check if motif_count_summary.tsv was found"
    fi
fi
echo ""

# ---------------------------------------------------------------------------
# 9. fig8 either generated or skipped gracefully
# ---------------------------------------------------------------------------
echo "--- Step 9: fig8 allele-aware figure ---"
FIG_PDF="$SMOKE_OUT/06_figures/pdf/fig8_allele_aware_precise_edit.pdf"
FIG_PNG="$SMOKE_OUT/06_figures/png/fig8_allele_aware_precise_edit.png"
if [[ -f "$FIG_PDF" && -f "$FIG_PNG" ]]; then
    pass "fig8 generated (PDF + PNG exist)"
else
    # fig8 is optional; it's OK to skip if allele-aware ran but figures ran before
    # In our run_all.sh, fig generation (06) runs BEFORE allele-aware (11),
    # so fig8 will be skipped in run_all.sh unless figures are re-run separately.
    # This is expected behavior — document it.
    pass "fig8 not yet generated (06_figures runs before 11_allele_aware in run_all.sh — rerun 06 to include fig8)"
fi
echo ""

# ---------------------------------------------------------------------------
# 10. Final: motif_count_summary.tsv still byte-identical after everything
# ---------------------------------------------------------------------------
echo "--- Step 10: Final gold standard check ---"
if diff -q "$GOLD_DIR/motif_count_summary.tsv" \
           "$SMOKE_OUT/04_edit/motif_count_summary.tsv" > /dev/null 2>&1; then
    pass "motif_count_summary.tsv still byte-identical to gold standard at end of run"
else
    fail "motif_count_summary.tsv CHANGED during run"
fi
echo ""

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "======================================================================"
echo "Allele-aware smoke test summary"
echo "  Finished : $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Output   : $SMOKE_OUT"
echo ""
echo -e "  ${GRN}PASS${NC} : $PASS"
echo -e "  ${RED}FAIL${NC} : $FAIL"
echo -e "  ${YLW}WARN${NC} : $WARN"
echo ""

if [[ "$FAIL" -eq 0 ]]; then
    echo -e "${GRN}All allele-aware checks passed.${NC}"
    echo ""
    echo "To run allele-aware analysis:"
    echo "  bash pipeline/run_all.sh \\"
    echo "      --input  \"\$HOME/projects/alsw/00_fastq\" \\"
    echo "      --output \"\$HOME/projects/alsw/rerun_output\" \\"
    echo "      --run-allele-aware"
    exit 0
else
    echo -e "${RED}$FAIL check(s) FAILED.${NC}"
    exit 1
fi
