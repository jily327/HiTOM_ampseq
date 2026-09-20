#!/usr/bin/env bash
# =============================================================================
# smoke_test_synthetic.sh  --  End-to-end check that needs no sequencing data.
#
# smoke_test.sh, smoke_test_alignment.sh and smoke_test_allele_aware.sh all
# read the real run from ~/projects/alsw/00_fastq, so they cannot run on a
# machine that does not hold it.  This test generates its own multiplexed
# FASTQ pair with a known composition (tests/make_synthetic_fastq.py), runs
# the whole pipeline including the optional modules, and checks the numbers
# the pipeline should produce for that input.
#
# It verifies behaviour, not the gold standard: reproducing
# outputs_example/motif_count_summary.tsv still requires the real FASTQs.
#
# Usage:
#   bash pipeline/tests/smoke_test_synthetic.sh [output_dir]
#
# Exit codes:
#   0 = all checks passed
#   1 = one or more checks failed
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TESTS_DIR="$SCRIPT_DIR/tests"
WORK_DIR="${1:-$(mktemp -d "${TMPDIR:-/tmp}/alsw_synth_XXXXXX")}"
FASTQ_DIR="$WORK_DIR/00_fastq"
OUT_DIR="$WORK_DIR/output"

RED='\033[0;31m'
GRN='\033[0;32m'
NC='\033[0m'

PASS=0
FAIL=0

pass() { echo -e "${GRN}PASS${NC}  $1"; ((PASS++)) || true; }
fail() { echo -e "${RED}FAIL${NC}  $1"; ((FAIL++)) || true; }

# check_field <file> <description> <python expression over row dicts>
check_value() {
    local label="$1" expected="$2" actual="$3"
    if [[ "$actual" == "$expected" ]]; then
        pass "$label = $actual"
    else
        fail "$label = $actual (expected $expected)"
    fi
}

echo "======================================================================"
echo "ALSW Portable Pipeline  —  Synthetic Smoke Test"
echo "  Pipeline dir : $SCRIPT_DIR"
echo "  Work dir     : $WORK_DIR"
echo "  Started      : $(date '+%Y-%m-%d %H:%M:%S')"
echo "======================================================================"
echo ""

# ---------------------------------------------------------------------------
# Interpreter
# ---------------------------------------------------------------------------
if [[ -z "${PYTHON:-}" ]]; then
    for _cand in python3 python python3.exe python.exe; do
        _path=$(command -v "$_cand" 2>/dev/null) || continue
        if "$_path" -c 'import sys; assert sys.version_info[:2] >= (3, 8)' \
                >/dev/null 2>&1; then
            PYTHON="$_path"
            break
        fi
    done
fi
if [[ -z "${PYTHON:-}" ]]; then
    fail "no working Python 3.8+ interpreter found"
    exit 1
fi
pass "Python interpreter: $($PYTHON --version 2>&1)"

# ---------------------------------------------------------------------------
# 1. Generate input
# ---------------------------------------------------------------------------
echo ""
echo "--- Generating synthetic FASTQs ---"
if "$PYTHON" "$TESTS_DIR/make_synthetic_fastq.py" "$FASTQ_DIR"; then
    pass "synthetic FASTQ pair written"
else
    fail "could not generate synthetic FASTQs"
    exit 1
fi

# ---------------------------------------------------------------------------
# 2. Run the pipeline
# ---------------------------------------------------------------------------
echo ""
echo "--- Running run_all.sh (with optional modules) ---"
if PYTHON="$PYTHON" bash "$SCRIPT_DIR/run_all.sh" \
        --input  "$FASTQ_DIR" \
        --output "$OUT_DIR" \
        --run-alignment \
        --run-allele-aware > "$WORK_DIR/run_all.log" 2>&1; then
    pass "run_all.sh exited 0"
else
    fail "run_all.sh returned non-zero exit code"
    tail -30 "$WORK_DIR/run_all.log"
    exit 1
fi

# ---------------------------------------------------------------------------
# 3. Expected outputs exist
# ---------------------------------------------------------------------------
echo ""
echo "--- Output files ---"
for rel in \
    "01_demux/demux_counts.tsv" \
    "02_qc/primer_dimer_qc_summary.csv" \
    "03_allele/allele_counts.tsv" \
    "04_edit/motif_count_summary.tsv" \
    "05_summary/final_sample_summary.csv" \
    "07_alignment/alignment_sample_summary.csv" \
    "10_comparison/motif_vs_alignment_comparison.csv" \
    "11_allele_aware/allele_aware_precise_edit_summary.csv"
do
    if [[ -f "$OUT_DIR/$rel" ]]; then
        pass "exists: $rel"
    else
        fail "missing: $rel"
    fi
done

# ---------------------------------------------------------------------------
# 4. Known-truth numbers
# ---------------------------------------------------------------------------
echo ""
echo "--- Values ---"

read_field() {
    # read_field <file> <delimiter> <match_col=value[,...]> <field>
    "$PYTHON" - "$@" <<'PY'
import csv, sys
path, delim, match, field = sys.argv[1:5]
pairs = [kv.split("=", 1) for kv in match.split(",")] if match else []
with open(path, newline="", encoding="utf-8") as fh:
    for row in csv.DictReader(fh, delimiter=("\t" if delim == "tsv" else ",")):
        if all(row.get(k) == v for k, v in pairs):
            print(row.get(field, ""))
            break
    else:
        print("<no matching row>")
PY
}

DEMUX="$OUT_DIR/01_demux/demux_counts.tsv"
check_value "demux assigned read pairs" "200" \
    "$(read_field "$DEMUX" tsv "sample=ePPEmax_ALSP_3dpi_pool" read_pairs)"
check_value "demux undetermined read pairs" "25" \
    "$(read_field "$DEMUX" tsv "sample=undetermined" read_pairs)"

MOTIF="$OUT_DIR/04_edit/motif_count_summary.tsv"
check_value "ALSP_chr01 desired %" "40.0000" \
    "$(read_field "$MOTIF" tsv "target=ALSP_chr01" desired_percent_among_motif_hits)"
check_value "ALSP_chr11 desired %" "70.0000" \
    "$(read_field "$MOTIF" tsv "target=ALSP_chr11" desired_percent_among_motif_hits)"
check_value "ALSP_core  desired %" "55.0000" \
    "$(read_field "$MOTIF" tsv "target=ALSP_core" desired_percent_among_motif_hits)"
check_value "ALSP_chr01 reads matching both motifs" "0" \
    "$(read_field "$MOTIF" tsv "target=ALSP_chr01" both_hits)"

# 02 must not report the undetermined bin
QC="$OUT_DIR/02_qc/primer_dimer_qc_summary.csv"
if grep -q "^undetermined," "$QC"; then
    fail "02_qc includes the undetermined bin"
else
    pass "02_qc excludes the undetermined bin"
fi

# 07: reads from the other homeolog must be reported as other_allele, and the
# desired percentage must still agree with the motif count.
ALIGN="$OUT_DIR/07_alignment/alignment_sample_summary.csv"
check_value "07 ALSP_chr01 other_allele reads" "100" \
    "$(read_field "$ALIGN" csv "target=ALSP_chr01" other_allele)"
check_value "07 ALSP_chr01 imprecise_PE reads" "0" \
    "$(read_field "$ALIGN" csv "target=ALSP_chr01" imprecise_PE)"
check_value "07 ALSP_chr01 desired %" "40.0000" \
    "$(read_field "$ALIGN" csv "target=ALSP_chr01" alignment_precise_desired_pct)"
check_value "07 ALSP_chr01 read pairs in sample" "200" \
    "$(read_field "$ALIGN" csv "target=ALSP_chr01" total_read_pairs_in_sample)"

# 10: motif and alignment must agree exactly on this clean input
CMP="$OUT_DIR/10_comparison/motif_vs_alignment_comparison.csv"
for tgt in ALSP_core ALSP_chr01 ALSP_chr11; do
    check_value "10 $tgt motif-vs-alignment difference" "0.0000" \
        "$(read_field "$CMP" csv "target=$tgt" difference_pct)"
done

# 11: allele-aware calls must match the motif counts, with nothing ambiguous
AA="$OUT_DIR/11_allele_aware/allele_aware_precise_edit_summary.csv"
check_value "11 chr01 desired %" "40.0000" \
    "$(read_field "$AA" csv "allele=chr01" precise_desired_pct_of_allele_assigned)"
check_value "11 chr11 desired %" "70.0000" \
    "$(read_field "$AA" csv "allele=chr11" precise_desired_pct_of_allele_assigned)"
check_value "11 ambiguous reads" "0" \
    "$(read_field "$AA" csv "allele=ambiguous" allele_assigned_reads)"

# ---------------------------------------------------------------------------
# 5. Guard rails
# ---------------------------------------------------------------------------
echo ""
echo "--- Guard rails ---"

# An unsynchronised FASTQ pair must fail loudly, not be silently truncated.
BAD_DIR="$WORK_DIR/00_fastq_truncated"
"$PYTHON" - "$FASTQ_DIR" "$BAD_DIR" <<'PY'
import gzip, shutil, sys
from pathlib import Path
src, dst = Path(sys.argv[1]), Path(sys.argv[2])
dst.mkdir(parents=True, exist_ok=True)
shutil.copy(src / "SYNTH_S1_L001_R1_001.fastq.gz", dst)
with gzip.open(src / "SYNTH_S1_L001_R2_001.fastq.gz", "rt", encoding="utf-8") as fh:
    lines = fh.readlines()
with gzip.open(dst / "SYNTH_S1_L001_R2_001.fastq.gz", "wt", encoding="utf-8") as fh:
    fh.writelines(lines[: 4 * 100])      # drop the tail of R2
PY
if "$PYTHON" "$SCRIPT_DIR/scripts/01_demux_and_qc.py" \
        --input "$BAD_DIR" --output "$WORK_DIR/truncated_out" \
        > "$WORK_DIR/truncated.log" 2>&1; then
    fail "01 accepted an unsynchronised R1/R2 pair"
elif grep -q "R1/R2 read counts differ" "$WORK_DIR/truncated.log"; then
    pass "01 rejects an unsynchronised R1/R2 pair"
else
    fail "01 failed on the truncated pair, but not with the expected message"
fi

# Config files holding non-ASCII text must load on any locale.
NONASCII="$WORK_DIR/targets_nonascii.json"
"$PYTHON" - "$SCRIPT_DIR/config/targets.json" "$NONASCII" <<'PY'
import json, sys
from pathlib import Path
cfg = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
cfg["_encoding_probe"] = "C→T edit — non-ASCII probe"
Path(sys.argv[2]).write_text(json.dumps(cfg, ensure_ascii=False, indent=2),
                             encoding="utf-8")
PY
if "$PYTHON" "$SCRIPT_DIR/scripts/04_edit_classification.py" \
        --input "$OUT_DIR" --output "$WORK_DIR/nonascii_out" \
        --config "$NONASCII" > "$WORK_DIR/nonascii.log" 2>&1; then
    pass "04 reads a config containing non-ASCII characters"
else
    fail "04 could not read a config containing non-ASCII characters"
    tail -5 "$WORK_DIR/nonascii.log"
fi

# ---------------------------------------------------------------------------
echo ""
echo "======================================================================"
echo "  Passed : $PASS"
echo "  Failed : $FAIL"
echo "  Work dir (kept for inspection): $WORK_DIR"
echo "======================================================================"

[[ $FAIL -eq 0 ]] || exit 1
exit 0
