#!/usr/bin/env bash
# =============================================================================
# smoke_test_alignment.sh  --  Smoke test for optional alignment modules.
#
# Tests:
#   1. Validated motif-count pipeline still passes (motif_count_summary.tsv
#      remains byte-identical to gold standard)
#   2. Alignment module (07) runs and produces expected outputs
#   3. CRISPResso2 prepare (08) runs without error (no CRISPResso2 required)
#   4. Comparison module (10) runs and produces expected outputs
#   5. CRISPResso2 is NOT required unless --run-crispresso2 is passed
#   6. make_new_hitom_project.py dry run
#
# Usage:
#   bash pipeline/tests/smoke_test_alignment.sh
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GOLD_DIR="$SCRIPT_DIR/outputs_example"
INPUT_DIR="$HOME/projects/alsw/00_fastq"
SMOKE_OUT="$HOME/projects/alsw/smoke_alignment_$(date +%Y%m%d_%H%M%S)"
PYTHON=$(command -v python3 || command -v python)

RED='\033[0;31m'; GRN='\033[0;32m'; YLW='\033[0;33m'; NC='\033[0m'
PASS=0; FAIL=0; WARN=0

pass() { echo -e "${GRN}PASS${NC}  $1"; ((PASS++)) || true; }
fail() { echo -e "${RED}FAIL${NC}  $1"; ((FAIL++)) || true; }
warn() { echo -e "${YLW}WARN${NC}  $1"; ((WARN++)) || true; }

echo "======================================================================"
echo "ALSW Alignment Module Smoke Test"
echo "  Pipeline    : $SCRIPT_DIR"
echo "  Input FASTQs: $INPUT_DIR"
echo "  Output      : $SMOKE_OUT"
echo "  Started     : $(date '+%Y-%m-%d %H:%M:%S')"
echo "======================================================================"
echo ""

# ---------------------------------------------------------------------------
# 1. Run the FULL validated pipeline (steps 01-06)
# ---------------------------------------------------------------------------
echo "--- Step 1: Run validated pipeline ---"
mkdir -p "$SMOKE_OUT"

if bash "$SCRIPT_DIR/run_all.sh" \
        --input  "$INPUT_DIR" \
        --output "$SMOKE_OUT" \
   > "$SMOKE_OUT/pipeline_run.log" 2>&1; then
    pass "Validated pipeline (01-06) completed cleanly"
else
    fail "Validated pipeline returned non-zero exit code"
    echo "  See: $SMOKE_OUT/pipeline_run.log"
    exit 1
fi

# Verify gold standard still matches
if diff -q "$GOLD_DIR/motif_count_summary.tsv" \
           "$SMOKE_OUT/04_edit/motif_count_summary.tsv" > /dev/null 2>&1; then
    pass "motif_count_summary.tsv byte-identical to gold standard"
else
    fail "motif_count_summary.tsv CHANGED after adding alignment module code"
fi
echo ""

# ---------------------------------------------------------------------------
# 2. Run alignment classifier (07) — ALSW target only for speed
# ---------------------------------------------------------------------------
echo "--- Step 2: Alignment classifier (07) ---"

if "$PYTHON" "$SCRIPT_DIR/scripts/07_alignment_edit_classifier.py" \
        --input  "$SMOKE_OUT" \
        --output "$SMOKE_OUT" \
        --target ALSW \
        --max-reads 500 \
   >> "$SMOKE_OUT/pipeline_run.log" 2>&1; then
    pass "07_alignment_edit_classifier.py ran without error (ALSW, 500 reads)"
else
    fail "07_alignment_edit_classifier.py returned non-zero exit code"
fi

EXPECTED_ALIGN=(
    "07_alignment/alignment_read_classifications.tsv"
    "07_alignment/alignment_sample_summary.csv"
    "07_alignment/alignment_indel_summary.csv"
    "07_alignment/alignment_vs_motif_comparison.csv"
)
for f in "${EXPECTED_ALIGN[@]}"; do
    if [[ -f "$SMOKE_OUT/$f" ]]; then
        pass "Exists: $f"
    else
        fail "Missing: $f"
    fi
done

# Verify motif_count_summary.tsv is unchanged (alignment must NOT touch it)
if diff -q "$GOLD_DIR/motif_count_summary.tsv" \
           "$SMOKE_OUT/04_edit/motif_count_summary.tsv" > /dev/null 2>&1; then
    pass "motif_count_summary.tsv unchanged after alignment step"
else
    fail "motif_count_summary.tsv was MODIFIED by alignment step — this is a bug"
fi
echo ""

# ---------------------------------------------------------------------------
# 3. CRISPResso2 input preparation (08) — no CRISPResso2 required
# ---------------------------------------------------------------------------
echo "--- Step 3: CRISPResso2 input prep (08) ---"

if "$PYTHON" "$SCRIPT_DIR/scripts/08_prepare_crispresso2_inputs.py" \
        --input  "$SMOKE_OUT" \
        --output "$SMOKE_OUT" \
   >> "$SMOKE_OUT/pipeline_run.log" 2>&1; then
    pass "08_prepare_crispresso2_inputs.py ran without error"
else
    fail "08_prepare_crispresso2_inputs.py returned non-zero exit code"
fi

if [[ -f "$SMOKE_OUT/08_crispresso2_inputs/CRISPRessoBatch_input.csv" ]]; then
    pass "CRISPRessoBatch_input.csv exists"
    # Verify it has FILL_IN placeholder (not a bug, just expected)
    if grep -q "FILL_IN_AMPLICON_SEQUENCE" \
            "$SMOKE_OUT/08_crispresso2_inputs/CRISPRessoBatch_input.csv"; then
        pass "Batch CSV has FILL_IN_AMPLICON_SEQUENCE placeholders (expected — no full amplicons provided)"
    else
        warn "No FILL_IN placeholder in batch CSV — unexpected for this dataset"
    fi
else
    fail "CRISPRessoBatch_input.csv not found"
fi

# Confirm CRISPResso2 is NOT required (script exits 0 when not installed)
CRISPRESSO_IN_PATH=$(command -v CRISPRessoBatch 2>/dev/null || echo "")
if [[ -z "$CRISPRESSO_IN_PATH" ]]; then
    if "$PYTHON" "$SCRIPT_DIR/scripts/09_run_crispresso2_batch.py" \
            --input  "$SMOKE_OUT" \
            --output "$SMOKE_OUT" \
       >> "$SMOKE_OUT/pipeline_run.log" 2>&1; then
        pass "09_run_crispresso2_batch.py exits 0 gracefully when CRISPResso2 not installed"
    else
        fail "09_run_crispresso2_batch.py should exit 0 when CRISPResso2 absent"
    fi
else
    warn "CRISPRessoBatch found in PATH ($CRISPRESSO_IN_PATH) — skipping not-installed test"
fi
echo ""

# ---------------------------------------------------------------------------
# 4. Motif vs alignment comparison (10)
# ---------------------------------------------------------------------------
echo "--- Step 4: Motif vs alignment comparison (10) ---"

if "$PYTHON" "$SCRIPT_DIR/scripts/10_compare_motif_vs_alignment.py" \
        --input  "$SMOKE_OUT" \
        --output "$SMOKE_OUT" \
   >> "$SMOKE_OUT/pipeline_run.log" 2>&1; then
    pass "10_compare_motif_vs_alignment.py ran without error"
else
    fail "10_compare_motif_vs_alignment.py returned non-zero exit code"
fi

for f in "10_comparison/motif_vs_alignment_comparison.csv" \
         "10_comparison/discrepancy_report.md"; do
    if [[ -f "$SMOKE_OUT/$f" ]]; then
        pass "Exists: $f"
    else
        fail "Missing: $f"
    fi
done
echo ""

# ---------------------------------------------------------------------------
# 5. make_new_hitom_project.py dry run
# ---------------------------------------------------------------------------
echo "--- Step 5: make_new_hitom_project.py ---"
NEW_PROJ="$SMOKE_OUT/test_new_project"
TEST_SAMPLES="$SMOKE_OUT/test_samples.csv"

# Create a minimal test sample sheet
cat > "$TEST_SAMPLES" << 'SAMPLEEOF'
sample_id,editor,target,dpi,replicate,forward_barcode,reverse_barcode
Test_TGT_3dpi,Editor1,TGT,3,rep1,F1,R-A
Test_TGT_6dpi,Editor1,TGT,6,rep1,F2,R-A
WT_TGT_3dpi,WT,TGT,3,rep1,F3,R-B
SAMPLEEOF

if "$PYTHON" "$SCRIPT_DIR/scripts/make_new_hitom_project.py" \
        --samples    "$TEST_SAMPLES" \
        --project-dir "$NEW_PROJ" \
        --project-name "Smoke Test Project" \
   >> "$SMOKE_OUT/pipeline_run.log" 2>&1; then
    pass "make_new_hitom_project.py ran without error"
else
    fail "make_new_hitom_project.py returned non-zero exit code"
fi

for f in "pipeline/config/sample_metadata.csv" \
         "pipeline/config/targets.json" \
         "pipeline/config/hitom_protocol_constants.json" \
         "pipeline/config/barcode_map.csv" \
         "pipeline/scripts/utils.py" \
         "pipeline/run_all.sh" \
         "PROJECT_README.md"; do
    if [[ -f "$NEW_PROJ/$f" ]]; then
        pass "New project: $f exists"
    else
        fail "New project: $f missing"
    fi
done

# Confirm targets.json template has FILL_IN
if grep -q "FILL_IN_WT_MOTIF_SEQUENCE" "$NEW_PROJ/pipeline/config/targets.json"; then
    pass "New project targets.json has FILL_IN placeholder (correct)"
fi
echo ""

# ---------------------------------------------------------------------------
# 6. Final summary
# ---------------------------------------------------------------------------
echo "======================================================================"
echo "Alignment smoke test summary"
echo "  Finished : $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Output   : $SMOKE_OUT"
echo ""
echo -e "  ${GRN}PASS${NC} : $PASS"
echo -e "  ${RED}FAIL${NC} : $FAIL"
echo -e "  ${YLW}WARN${NC} : $WARN"
echo ""

if [[ "$FAIL" -eq 0 ]]; then
    echo -e "${GRN}All alignment module checks passed.${NC}"
    exit 0
else
    echo -e "${RED}$FAIL check(s) FAILED.${NC}"
    exit 1
fi
