# CRISPResso2 Local Usage Guide

## What CRISPResso2 does (vs this pipeline)

| | Motif-count pipeline (validated) | CRISPResso2 |
|---|---|---|
| Method | Exact substring search | Read alignment to amplicon |
| Indel detection | No | Yes |
| Imprecise PE detection | No | Yes |
| Requires reference amplicon | No | Yes |
| Validated for this dataset | Yes | Not yet |
| Speed | Fast (pure Python) | Moderate |
| Output | TSV counts | HTML reports + allele tables |

The motif-count pipeline is the **validated gold standard** for this project.
CRISPResso2 provides complementary alignment-based analysis.

## Installing CRISPResso2

### Option A — conda (recommended)

```bash
conda create -n crispresso2 -c conda-forge -c bioconda crispresso2 -y
conda activate crispresso2
CRISPRessoBatch --version
```

### Option B — pip

```bash
pip install CRISPResso2
```

### Option C — Docker

```bash
docker pull pinellolab/crispresso2
alias CRISPRessoBatch="docker run -v $(pwd):/DATA -w /DATA pinellolab/crispresso2 CRISPRessoBatch"
```

## Preparing inputs

```bash
python3 pipeline/scripts/08_prepare_crispresso2_inputs.py \
    --input  "$HOME/projects/alsw/rerun_output" \
    --output "$HOME/projects/alsw/rerun_output"
```

This creates `08_crispresso2_inputs/CRISPRessoBatch_input.csv`.

**Before running CRISPResso2**, edit the CSV and replace
`FILL_IN_AMPLICON_SEQUENCE` with the correct full amplicon sequence for
each target.

## Running CRISPResso2

### Via pipeline wrapper

```bash
bash pipeline/run_all.sh \
    --input  "$HOME/projects/alsw/00_fastq" \
    --output "$HOME/projects/alsw/rerun_output" \
    --run-crispresso2
```

### Manually

```bash
conda activate crispresso2
CRISPRessoBatch \
    --batch_settings "$HOME/projects/alsw/rerun_output/08_crispresso2_inputs/CRISPRessoBatch_input.csv" \
    --output_folder  "$HOME/projects/alsw/rerun_output/09_crispresso2_results" \
    --n_processes 4 \
    --min_frequency_alleles_around_cut_to_plot 0.05
```

## Interpreting CRISPResso2 outputs

Key output files per sample:

| File | Contents |
|------|---------|
| `CRISPResso_on_*/Alleles_frequency_table.zip` | All observed allele sequences + frequencies |
| `CRISPResso_on_*/CRISPResso2_quantification_of_editing_frequency.txt` | Summary: % NHEJ, % HDR, % mixed |
| `CRISPResso_on_*/CRISPRessoBatch_alleles_frequency_table.txt` | Batch-level allele summary |

Key metrics:
- **% HDR** (homology-directed repair) ≈ precise desired edit fraction
- **% NHEJ** ≈ indel fraction
- **% Unmodified** ≈ WT fraction

### Relationship to motif-count results

| Metric | Motif-count | CRISPResso2 |
|--------|-------------|-------------|
| "Editing %" | desired / (WT + desired) motif hits | HDR / (total mapped reads) |
| Denominator | informative reads (WT + desired) | all mapped reads |
| Indels | not counted | counted as NHEJ |

The denominators differ, so percentages will not be identical. Use
`10_compare_motif_vs_alignment.py` to compare the two estimates.

## What is NOT shown by CRISPResso2 in this dataset

- Full amplicon sequences are required but not yet provided for ALSW, ALSP,
  and EPSPS targets. Fill in `reference_amplicon_sequence` in `targets.json`
  and re-run `08_prepare_crispresso2_inputs.py`.
- Without the correct amplicon sequence, CRISPResso2 will fail or produce
  incorrect results.

## Validation note

CRISPResso2 outputs for this dataset have **not been cross-validated** against
the motif-count pipeline. Always compare both results using
`10_compare_motif_vs_alignment.py` before drawing conclusions.
