# ALSW HiTOM Amplicon-seq Pipeline

Portable, reproducible pipeline for prime-editing efficiency analysis from
HiTOM multiplexed amplicon sequencing data.

---

## What this pipeline does

This pipeline takes raw paired-end FASTQ files from a HiTOM-barcoded amplicon
sequencing experiment, demultiplexes them by dual barcode, counts WT and
desired-edit motifs in each sample, and produces final summary tables and
figures.

The pipeline is a direct port of a validated ad-hoc analysis run on 2026-06-01.
The gold-standard output (`outputs_example/motif_count_summary.tsv`) was
produced by that session and is used to verify reproducibility.

---

## Quick start

```bash
# Copy raw FASTQs to WSL-local directory (one-time step)
cp /mnt/d/... ~/projects/alsw/00_fastq/

# Run full pipeline
bash pipeline/run_all.sh \
    --input  "$HOME/projects/alsw/00_fastq" \
    --output "$HOME/projects/alsw/rerun_output"
```

For figures, first install dependencies:
```bash
pip install matplotlib pandas numpy
```

`run_all.sh` picks an interpreter by running each candidate, so a shim that
is on `PATH` but does not work (such as the Microsoft Store `python3` stub on
Windows) is skipped. Set `PYTHON=/path/to/python` to choose one yourself.

Steps 01 to 05 and 07 to 11 need only the Python standard library. Steps 06
and the scatter plot in step 10 need matplotlib, pandas and numpy; without
them those steps are skipped and the rest of the pipeline still completes.

---

## Input files

### Required
| File | Description |
|------|-------------|
| `*_R1_*.fastq.gz` | Raw paired-end R1 reads (one file per run) |
| `*_R2_*.fastq.gz` | Raw paired-end R2 reads (one file per run) |

### Config (provided in `pipeline/config/`, do not need to change for same experiment)
| File | Description |
|------|-------------|
| `sample_metadata.csv` | Sample-to-barcode assignments for this run |
| `barcode_map.csv` | Resolved barcode sequences (derived from metadata) |
| `targets.json` | WT and desired motif sequences for each target |
| `hitom_protocol_constants.json` | Fixed HiTOM library design constants |

---

## Expected FASTQ naming

The pipeline expects **exactly one R1 and one R2 file** in the input directory,
named following standard Illumina conventions:
```
*_R1_001.fastq.gz   or   *_R1.fastq.gz
*_R2_001.fastq.gz   or   *_R2.fastq.gz
```

---

## Where default config files live

All default configs are in `pipeline/config/`.  The pipeline finds them
automatically relative to `run_all.sh`.  You can override any config with
`--metadata`, `--config`, or `--protocol-constants` flags.

---

## HiTOM protocol constants (`hitom_protocol_constants.json`)

These are **fixed by the HiTOM library design** and should not change between
experiments unless the library design changes.

### Bridge sequences
| Name | Sequence |
|------|----------|
| Forward bridge | `ggagtgagtacggtgtgc` |
| Reverse bridge | `gagttggatgctggatgg` |

These bridges appear as tail sequences when the amplicon is shorter than the
read length:
- **R1 tail bridge** = RC(reverse_bridge) = `CCATCCAGCATCCAACTC`
- **R2 tail bridge** = RC(forward_bridge) = `GCACACCGTACTCACTCC`

### Standard Qi-lab forward barcodes (12-plex)
| Name | Seq | Name | Seq | Name | Seq |
|------|-----|------|-----|------|-----|
| F1 | GCGT | F5 | GCTC | F9  | ATAC |
| F2 | GTAG | F6 | AGTC | F10 | CACA |
| F3 | ACGC | F7 | CGAC | F11 | GTGC |
| F4 | CTCG | F8 | GATG | F12 | ACTA |

### Standard Qi-lab reverse barcodes (8-plex)
| Name | Seq | Name | Seq |
|------|-----|------|-----|
| R-A | GCGT | R-E | GCTC |
| R-B | GTAG | R-F | AGTC |
| R-C | ACGC | R-G | CGAC |
| R-D | CTCG | R-H | GATG |

---

## How sample-specific barcode assignment works

`pipeline/config/sample_metadata.csv` maps each sample to a barcode pair by
**name** (e.g. `F1`, `R-E`).  The pipeline looks up the actual 4-bp sequences
from `hitom_protocol_constants.json` at runtime.

This means you can use barcode names instead of raw sequences, and the pipeline
will validate them against the known dictionary and print an error if an unknown
barcode name is used.

You can also use raw 4-bp sequences directly in the metadata (e.g. `GCGT`
instead of `F1`); the pipeline accepts both.

---

## How barcode demultiplexing works

Reads are assigned by exact match of:
- R1 positions 4–8 (forward barcode)
- R2 positions 4–8 (reverse barcode)

Reads with no matching barcode pair → `undetermined`.

This is exact-match only.  No mismatch tolerance is applied (consistent with
the validated analysis).

---

## How 27-bp leader trimming works

After barcode assignment, the first 27 bases of both R1 and R2 are removed.

The 27-bp leader structure:
```
[4 random bp] + [4 bp barcode] + [~19 bp bridge/primer]
```

This leaves only the amplicon sequence for downstream motif counting.

---

## How tail bridge clipping works

When the amplicon is shorter than the read length, R1 reads through the amplicon
and into the reverse bridge sequence, and vice versa for R2.

After leader trimming, any occurrence of the tail bridge sequence is clipped:
- R1 tail bridge: `CCATCCAGCATCCAACTC` (RC of reverse bridge)
- R2 tail bridge: `GCACACCGTACTCACTCC` (RC of forward bridge)

---

## How motif counting works

For each target × applicable sample:

1. Read R1 and R2 reads in parallel from `demux_clean/`
2. Build combined search string:
   ```
   combined = R1 + N + R2 + N + RC(R1) + N + RC(R2)
   ```
3. Check for exact substring match of WT motif and desired motif
   (including reverse complement of each motif in the combined string)
4. Count:
   - `wt_motif_hits`: reads with WT motif
   - `desired_motif_hits`: reads with desired motif
   - `both_hits`: reads with both
   - `neither_hits`: reads with neither
   - `informative_reads` = wt_hits + desired_hits
   - `desired_percent_among_motif_hits` = desired_hits / informative × 100

This is the **validated denominator**.

---

## What is validated

The following outputs are validated against the gold-standard
`outputs_example/motif_count_summary.tsv`:

- `04_edit/motif_count_summary.tsv` — must match exactly (to 4 decimal places)
- `01_demux/demux_counts.tsv` — read counts per sample

The smoke test (`tests/smoke_test.sh`) verifies these automatically. It needs the
real run in `~/projects/alsw/00_fastq`.

### Running without the sequencing data

`tests/smoke_test_synthetic.sh` needs no external files. It generates a small
multiplexed FASTQ pair of known composition, runs the whole pipeline including
the optional modules, and checks the numbers that input must produce:

```bash
bash pipeline/tests/smoke_test_synthetic.sh
```

Use it to check that the code works on a new machine or after a change. It
does not reproduce the gold standard; only `tests/smoke_test.sh` against the
real FASTQs does that.

---

## What is new QC/convenience code

| Script | Status |
|--------|--------|
| `01_demux_and_qc.py` | **VALIDATED** (port of demux_hitom.py + trim_tail_bridges.py) |
| `02_primer_dimer_filter.py` | **NEW** — informational QC only; does not filter reads |
| `03_allele_assignment.py` | **VALIDATED** for ALSP chr01/chr11; ALSW single-allele |
| `04_edit_classification.py` | **VALIDATED** — produces gold-standard motif_count_summary.tsv |
| `05_make_summary_tables.py` | **DERIVED** — aggregates validated tables |
| `06_make_figures.py` | **NEW** — convenience figures; no new biological conclusions |

---

## What is NOT validated

### Indel and imprecise-PE counting
The `edit_counts.tsv` table contains columns `indel_reads` and
`imprecise_PE_reads` filled with `NOT_VALIDATED_alignment_required`.

**Indel and imprecise-PE classification requires alignment-based calling**
(e.g. CRISPResso2 or BWA + variant calling pipeline).  The current validated
analysis is motif-based and binary: WT motif vs desired motif.

These stub columns exist to keep the table schema forward-compatible for when
alignment-based calling is implemented.

### ALSW chr01/chr11 allele split
The validated ALSW analysis uses a **single WT/desired motif pair** with no
chromosome-specific allele split.  If allele-discriminating flanking SNPs for
ALSW chr01 and chr11 are identified experimentally, they can be added to
`targets.json` under the `applicable_to: ["ALSW"]` group without any code
changes.

### Primer-dimer filtering
The `02_primer_dimer_filter.py` step computes short-read fractions but does
**not** remove short reads from the motif-counting path.  The validated motif
counts in `motif_count_summary.tsv` use all reads from `demux_clean/` (this is
identical to the validated original).

---

## Output descriptions

```
<output>/
├── 01_demux/
│   ├── demux_trimmed/        per-sample FASTQs after barcode demux + 27-bp trim
│   ├── demux_clean/          per-sample FASTQs after tail-bridge clip  ← used downstream
│   ├── demux_counts.tsv      read pairs per sample
│   ├── tail_trim_summary.tsv tail-bridge clip counts
│   ├── barcode_assignment_summary.csv  overall assignment rate
│   ├── read_length_summary.csv         post-trim length stats
│   ├── barcode_map_resolved.csv        barcode names resolved to sequences
│   └── md5_summary.txt                 MD5 check results
├── 02_qc/
│   ├── primer_dimer_qc_summary.csv     short-read fractions (informational)
│   ├── per_sample_length_histogram.csv R1 length histograms
│   └── potential_dimer_reads_summary.csv
├── 03_allele/
│   ├── allele_counts.tsv               per-target allele counts
│   └── allele_specific_summary.csv     ALSP chr01/chr11 split
├── 04_edit/
│   ├── motif_count_summary.tsv  ← GOLD STANDARD OUTPUT
│   └── edit_counts.tsv          extended table with NOT_VALIDATED stub columns
├── 05_summary/
│   ├── final_sample_summary.csv
│   ├── per_condition_summary.csv
│   ├── allele_specific_summary.csv
│   ├── primer_dimer_qc_summary.csv
│   ├── wt_background_summary.csv
│   └── sensitivity_denominators.csv
└── 06_figures/
    ├── 01_editing_efficiency_by_sample.{png,pdf}
    ├── 02_editor_comparison_ALSW.{png,pdf}
    ├── 03_dpi_heatmap.{png,pdf}
    ├── 04_allele_composition_ALSP.{png,pdf}
    ├── 05_primer_dimer_fraction.{png,pdf}
    ├── 06_wt_background.{png,pdf}
    └── 07_target_edit_composition.{png,pdf}
```

---

## Validated motif-count analysis

The **only validated analysis** in this pipeline is the motif-count path:
Steps 01 → 02 → 03 → 04 → 05 → 06. The gold-standard output is:

```
<output>/04_edit/motif_count_summary.tsv
```

All biological conclusions (editing %, WT background, allele-specific counts)
cited in reports or manuscripts should come from this file.

---

## Optional allele-aware precise edit quantification

Step 11 (`11_allele_aware_precise_edit.py`) assigns reads to individual
homeologs (e.g. ALSP chr01 vs chr11) using allele-discriminating SNPs,
then calls edit status per allele.

```bash
bash pipeline/run_all.sh \
    --input  "$HOME/projects/alsw/00_fastq" \
    --output "$HOME/projects/alsw/rerun_output" \
    --run-allele-aware
```

**ALSP is configured** — the allele-discriminating SNP at position 13 of the
29 bp context motif (chr01: T; chr11: C) is derived from the validated
ALSP_chr01 / ALSP_chr11 motif pair and stored in `config/targets.json`.

**ALSW is unconfigured** — no validated allele-discriminating SNPs are
available. The script skips ALSW gracefully with a warning.

**Validation status: NOT VALIDATED.** Results are exploratory. Compare
against `allele_aware_vs_motif_comparison.csv` and validate against controls
before reporting.

See `docs/allele_aware_precise_edit.md` for full methodology.

---

## Optional alignment-based indel/PE classification

Step 07 (`07_alignment_edit_classifier.py`) adds alignment-based read
classification using pure-Python Smith-Waterman.

```bash
bash pipeline/run_all.sh \
    --input  "$HOME/projects/alsw/00_fastq" \
    --output "$HOME/projects/alsw/rerun_output" \
    --run-alignment
```

**Validation status: NOT VALIDATED.**
Treat alignment outputs as exploratory until confirmed against controls.

By default step 07 classifies at most 3000 read pairs per sample and target,
drawn at random from the whole file (`--sampling random --seed 0`), while step
04 always uses every read. Each summary row records how many reads the sample
holds and how many were classified, so the two methods can be compared
honestly. Reads belonging to the other homeolog of an allele-specific target
are reported as `other_allele` rather than as an editing outcome.

Key limitation: full amplicon reference sequences are not currently provided
in `targets.json` (`reference_amplicon_sequence: null`). Alignment uses the
short WT/desired motif (~15–29 bp) as reference. Results describe the motif
window only.

See `docs/alignment_classifier_notes.md` for full details.

---

## Optional CRISPResso2 local analysis

Steps 08–09 prepare inputs for and optionally run CRISPResso2.

```bash
bash pipeline/run_all.sh \
    --input  "$HOME/projects/alsw/00_fastq" \
    --output "$HOME/projects/alsw/rerun_output" \
    --prepare-crispresso2   # prepares batch CSV; fills in PLACEHOLDER for missing amplicons
    # --run-crispresso2     # also runs CRISPResso2 if installed
```

CRISPResso2 is **not installed** in the pipeline by default. See
`docs/crispresso2_local_usage.md` for installation instructions.

**Important**: amplicon sequences (`reference_amplicon_sequence` in
`targets.json`) must be provided before CRISPResso2 can run.

---

## When to trust which output

| Output | Trust level | Use for |
|--------|-------------|--------|
| `04_edit/motif_count_summary.tsv` | ✅ **Validated** | All biological conclusions |
| `05_summary/final_sample_summary.csv` | ✅ Derived from validated | Summary tables |
| `06_figures/` | ✅ Derived from validated | Figures |
| `07_alignment/alignment_sample_summary.csv` | ⚠️ Exploratory | Method comparison only |
| `09_crispresso2_results/` | ⚠️ Third-party, not cross-validated | Independent verification |
| `10_comparison/motif_vs_alignment_comparison.csv` | ⚠️ Comparison only | Flagging discrepancies |

---

## Future run template generator

To set up the pipeline for a new HiTOM experiment:

```bash
# Create a sample sheet CSV with columns:
# sample_id, target, editor, dpi, replicate, forward_barcode, reverse_barcode
vim my_new_samples.csv

# Generate new project
python3 pipeline/scripts/make_new_hitom_project.py \
    --samples      my_new_samples.csv \
    --project-dir  ~/projects/new_experiment \
    --project-name "New HiTOM Experiment"

# Then fill in motif sequences in the generated targets.json
vim ~/projects/new_experiment/pipeline/config/targets.json

# Run
bash ~/projects/new_experiment/pipeline/examples/run_this_project.sh
```

Barcode columns accept either Hi-TOM names (`F1`, `R-A`) or raw 4-bp
sequences (`GCGT`, `GCTC`). The pipeline validates names against the
built-in dictionary.

---

## Figure generation

Figures are automatically regenerated by `run_all.sh` (Step 06). They require
`matplotlib`, `pandas`, and `numpy`; Steps 01–05 work without them.

```bash
pip install matplotlib pandas numpy
```

All figures are saved in three formats under `<output>/06_figures/`:

```
06_figures/
├── pdf/       vector PDFs  (for manuscripts)
├── svg/       vector SVGs  (for slides / editing)
├── png/       300 dpi raster PNGs  (for quick viewing)
├── source_data/  source CSV for each figure
└── figure_legends.md  auto-generated draft legends
```

You can also run figure generation alone:
```bash
python3 pipeline/scripts/06_make_figures.py \
    --input  "$HOME/projects/alsw/rerun_output" \
    --output "$HOME/projects/alsw/rerun_output" \
    --format pdf,png,svg \
    --dpi 300
```

### Which figures are validated result figures

| Figure | Contents | Validation status |
|--------|----------|------------------|
| Fig 1 | Editing efficiency by editor and target | **Validated** — sourced from motif_count_summary.tsv |
| Fig 2 | ePPEmax vs PE6c comparison (ALSW) | **Validated** |
| Fig 3 | 3 dpi vs 6 dpi effect (ALSW) | **Validated** |
| Fig 4 | WT background / false-positive | **Validated** |
| Fig 5 | Primer-dimer / short-read QC | **QC only** — informational, not validated analysis |
| Fig 6 | Motif composition stacked bar | **Validated** counts; visualisation only |
| Fig 7 | ALSP allele-specific (chr01/chr11) | **Validated** for ALSP; ALSW not shown (no allele split) |

### Indel and imprecise-PE are not shown in figures

Indel and imprecise prime-editing categories require alignment-based classification
(e.g. CRISPResso2). They are not plotted because they are not available from
this motif-count pipeline. Stub columns exist in `edit_counts.tsv` for
forward compatibility.

### ALSW is not allele-split in figures

Figure 7 only shows ALSP chr01/chr11. ALSW is treated as single-allele
because no validated allele-discriminating flanking motifs exist for that
target in this analysis. Add them to `targets.json` when available.

---

## Known limitations

1. **No indel or imprecise-PE calling** — alignment-based classification not
   implemented.  Stub columns in `edit_counts.tsv` are zero-filled placeholders.

2. **ALSW is single-allele** — No validated allele-discriminating motifs for
   ALSW chr01/chr11 exist in this analysis.

3. **Exact barcode match only** — Mismatch-tolerant demux is not used (consistent
   with the validated pipeline).

4. **Informative reads denominator** — The main editing percentage uses
   `informative_reads = wt_hits + desired_hits` as denominator, not total read
   pairs.  This means reads with neither WT nor desired motif (usually the
   majority, due to short/non-specific products) are excluded from the
   denominator.  This is the validated approach.

5. **Primer-dimer reads are not filtered** — Short reads are counted but not
   removed.  This is intentional to preserve the validated motif-count numbers.

---

## Regenerated vs static outputs

| Output | Regenerated from raw FASTQ | Static |
|--------|---------------------------|--------|
| `04_edit/motif_count_summary.tsv` | ✅ Yes | — |
| `01_demux/demux_counts.tsv` | ✅ Yes | — |
| `05_summary/*.csv` | ✅ Yes (derived) | — |
| `06_figures/*.png/pdf` | ✅ Yes | — |
| `outputs_example/motif_count_summary.tsv` | — | ✅ Gold standard |
