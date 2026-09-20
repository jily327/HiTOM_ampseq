# Allele-Aware Precise Edit Quantification

## Why allele-aware counting is needed

In polyploid organisms such as *Nicotiana benthamiana* (a near-allotetraploid)
and hexaploid wheat (*Triticum aestivum*), the target gene exists as homeologous
copies on different chromosome sets. For example, the *ALS* gene has homeologs
on chromosomes 1 and 11 in *N. benthamiana*.

A standard HiTOM amplicon experiment co-amplifies **both homeologs** from a
single primer pair. The motif-count pipeline (`04_edit_classification.py`)
counts reads matching the desired motif regardless of which homeolog they came
from. This is the validated, conservative estimate.

Allele-aware analysis asks the additional question:
> *Of the chr01 reads, what fraction carry the precise desired edit?
>  Of the chr11 reads, what fraction carry the precise desired edit?*

This matters because:
- Editing efficiency may differ between homeologs.
- A desired edit in one homeolog does not imply editing of the other.
- The two homeologs may differ in baseline susceptibility to the guide RNA.

## How allele assignment works

**Key principle: allele assignment uses ONLY positions that differ between
homeologs and are OUTSIDE the edit window.**

### Step 1 — Find the allele-discriminating anchor

The algorithm searches for a common anchor sequence present in both homeologs.
For ALSP:

```
chr01: TTGTTGCTATAAC T GGTCAAGTG C CACGT
                    ^            ^
                  pos 13       pos 23
               allele SNP     edit site
                  (T=chr01)   (C=WT; T=desired)

chr11: TTGTTGCTATAAC C GGTCAAGTG C CACGT
                    ^
                  (C=chr11)
```

The anchor `TTGTTGCTATAAC` (positions 0–12 of the 29 bp motif) is searched in
the combined read string (R1 + N + R2 + N + RC(R1) + N + RC(R2)).

### Step 2 — Check allele-discriminating SNP

The base immediately after the anchor (position 13) is the allele-discriminating
SNP:
- **T** → assign to **chr01**
- **C** → assign to **chr11**
- Other / not found → **low_information** or **ambiguous**

### Step 3 — Call edit status (assigned reads only)

The edit site is at anchor_end + 10 (position 23 in the 29 bp motif):
- **C** → **wt** (not edited)
- **T** → **precise_desired** (correctly edited)
- Other → **non_desired**
- Not covered → **ambiguous_edit**

## Why the edit base is excluded from allele assignment

The desired edit (C→T at position 23) changes the edit position but does **not**
distinguish chr01 from chr11 — both homeologs have C at position 23 in the
WT state. If we allowed the edit position to influence allele scoring, a
precisely edited read with T at position 23 could be misassigned.

The `exclude_from_allele_assignment: true` flag in `targets.json` ensures this.
An additional `edit_window` (± 3 bp around position 23) is also excluded as a
conservative safeguard.

## How ambiguous reads are handled

| Category | Meaning | Action |
|----------|---------|--------|
| `chr01` | SNP unambiguously matches chr01 | Counted in chr01 totals |
| `chr11` | SNP unambiguously matches chr11 | Counted in chr11 totals |
| `ambiguous` | SNP score tied; or SNP base unexpected | Excluded from allele-specific denominators |
| `low_information` | Anchor not found in read | Excluded; reported in QC table |

A read is called **ambiguous** if the allele margin (score_chr01 − score_chr11)
is less than `--min-allele-margin` (default: 1). For a single SNP, this means
the winning allele must have at least 1 more matching SNP than the other, which
effectively means a clear single-SNP match.

## Denominators

Two denominators are reported:

1. **`precise_desired_pct_of_allele_assigned`** — precise edits divided by all
   reads assigned to that allele (chr01 or chr11 reads, regardless of edit
   status). This is the most direct allele-specific estimate.

2. **`precise_desired_pct_of_total_informative`** — precise edits for that
   allele divided by all allele-assigned reads (both alleles combined). This
   is useful for comparing the per-allele contribution to total editing.

The **validated** editing percentage from the motif-count pipeline uses
`informative_reads = wt_hits + desired_hits` as denominator. These are
different denominators and numbers will not match.

## How to add a new allele-aware target to targets.json

Add an entry to the `allele_aware_targets` array:

```json
{
  "target_name":          "MY_TARGET",
  "analysis_type":        "allele_aware_precise_edit",
  "applicable_to":        ["MY_TARGET_GENE"],
  "status":               "configured",
  "allele_common_anchor": "COMMON_PREFIX_SEQUENCE",
  "alleles": {
    "allele_A": {
      "reference_sequence": null,
      "allele_snps": [
        {
          "snp_offset_from_anchor_end": 0,
          "expected_base": "A",
          "description": "allele_A-specific base at SNP position"
        }
      ]
    },
    "allele_B": {
      "reference_sequence": null,
      "allele_snps": [
        {
          "snp_offset_from_anchor_end": 0,
          "expected_base": "G",
          "description": "allele_B-specific base at SNP position"
        }
      ]
    }
  },
  "edit": {
    "name":              "my_edit",
    "offset_from_snp":   N,
    "wt_base":           "C",
    "desired_base":      "T",
    "exclude_from_allele_assignment": true
  },
  "min_informative_snps": 1,
  "notes":               "Description of target."
}
```

Then run:
```bash
bash pipeline/run_all.sh \
    --input  "$HOME/projects/alsw/00_fastq" \
    --output "$HOME/projects/my_output" \
    --run-allele-aware
```

## Limitations and validation requirements

1. **Single SNP for ALSP** — allele assignment is based on a single
   allele-discriminating position. A single sequencing error at that position
   could cause mis-assignment. Results with `--min-informative-snps 1` are
   less robust than multi-SNP discrimination.

2. **No full amplicon reference** — the algorithm searches for the anchor in
   reads but does not perform full-read alignment. It only works if the anchor
   sequence is present in the read. Primer-dimer and off-target short reads
   that do not contain the anchor are classified as `low_information`.

3. **Not validated** — these results have not been compared to known-mixture
   controls or orthogonal methods. Treat as exploratory.

4. **ALSW is unconfigured** — no allele-discriminating SNP positions are
   available for ALSW homeologs. ALSW is reported as `status: unconfigured`
   and skipped gracefully.

5. **Comparison with motif-count** — allele-aware percentages use different
   denominators than motif_count_summary.tsv and will not match numerically.
   Cross-check using `allele_aware_vs_motif_comparison.csv`.

## Output file descriptions

| File | Description |
|------|-------------|
| `allele_aware_precise_edit_summary.csv` | Per-sample × per-allele counts and percentages |
| `read_level_allele_edit_calls.tsv` | Per-read allele assignment and edit call |
| `allele_assignment_qc.csv` | Assignment rates per sample |
| `allele_aware_vs_motif_comparison.csv` | Side-by-side with validated motif counts |
