# Alignment-Based Edit Classifier — Notes

## What this module does

`07_alignment_edit_classifier.py` uses pure-Python Smith-Waterman local
alignment to classify reads beyond the binary WT/desired-motif categories
produced by the validated motif-count pipeline.

## Validation status

**NOT VALIDATED.** This module is newly implemented and has not been checked
against known reference standards or biological controls. All outputs should
be treated as exploratory until confirmed independently.

The validated gold standard for this project remains:

> `04_edit/motif_count_summary.tsv`

Do not use alignment-derived counts in publications or reports without
independent validation.

## Important limitation: no full amplicon reference

The current `targets.json` has `reference_amplicon_sequence: null` for all
targets because the full amplicon sequences were not provided during pipeline
setup. This means:

- Alignment is performed against the **short WT/desired motif sequences**
  (~15–29 bp), not the full amplicon.
- Indels and substitutions **outside the motif window are not detected.**
- Results describe the motif window only.

To enable full-amplicon alignment, add the complete amplicon sequence to
`targets.json`:

```json
"reference_amplicon_sequence": "ATCGATCG...full sequence...ATCGATCG"
```

## Classification categories

| Category | Meaning |
|----------|---------|
| `wt` | Exact match to WT motif (consistent with motif-count) |
| `precise_desired` | Exact match to desired motif, no additional indels |
| `desired_plus_indel` | Desired motif present AND indel detected in alignment |
| `imprecise_PE` | Partial desired change (some but not all edit positions changed) |
| `indel` | Indel detected in motif window, no exact motif match |
| `other_substitution` | Substitution not matching WT or desired at edit positions |
| `no_motif_align` | Alignment score below identity threshold |
| `short_read` | Read too short for alignment (< --min-read-len, default 30 bp) |

## Speed notes

- Smith-Waterman is O(n × m) per alignment, where n = read length (~100 bp),
  m = motif length (~15–29 bp).
- With `--max-reads 3000` (default), alignment runs in ~30–60 seconds per
  sample for 18 samples.
- Set `--max-reads 0` to process all reads (much slower).

## How to validate this output

1. Compare `alignment_vs_motif_comparison.csv` against `motif_count_summary.tsv`.
   The desired % estimates should broadly agree for the ALSW target.
2. Manually inspect reads classified as `indel` or `imprecise_PE` to confirm
   the alignment is biologically plausible.
3. Run on a synthetic control: a FASTQ containing known WT and edited reads
   at known ratios, and verify the alignment classifier recovers the correct
   proportions.

## Known caveats

- The short-read fraction (75–96% in this dataset) means most reads are
  classified as `short_read` and skipped. This is expected.
- The "neither" reads in the motif-count pipeline include both primer dimers
  (short, correctly excluded) and potentially real amplicon reads with indels
  (longer, classified by alignment).
- Imprecise PE detection is based on alignment to the short motif, so context
  outside the motif is not examined.
