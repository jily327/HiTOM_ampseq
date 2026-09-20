#!/usr/bin/env python3
"""
make_synthetic_fastq.py  --  Build a HiTOM FASTQ pair with known answers.

The other smoke tests need the real sequencing run in ~/projects/alsw/00_fastq,
so they cannot run on a machine that does not hold the data.  This generator
writes a small multiplexed FASTQ pair whose composition is fixed and known, so
that smoke_test_synthetic.sh can check the pipeline end to end anywhere.

Composition (sample ePPEmax_ALSP_3dpi_pool, barcodes F1 / R-E):

    ALSP chr01  WT       60 read pairs
    ALSP chr01  desired  40 read pairs
    ALSP chr11  WT       30 read pairs
    ALSP chr11  desired  70 read pairs
    unknown forward barcode  25 read pairs  -> undetermined

Expected results:

    demux                    200 assigned / 25 undetermined
    ALSP_chr01  desired      40.0000 %   (40 / (60 + 40))
    ALSP_chr11  desired      70.0000 %   (70 / (30 + 70))
    ALSP_core   desired      55.0000 %   (110 / 200, alleles pooled)

Every read is a clean, error-free molecule, so the motif counts are exact.
Reads are emitted in a fixed shuffled order (seeded) to keep runs comparable.

Usage:
    python tests/make_synthetic_fastq.py <output_dir>
"""

import gzip
import random
import sys
from pathlib import Path

SEED = 20260919

# Hi-TOM protocol constants (must match config/hitom_protocol_constants.json)
FORWARD_BARCODE = "GCGT"           # F1
REVERSE_BARCODE = "GCTC"           # R-E
UNKNOWN_BARCODE = "TTTT"           # not in the 12-plex set -> undetermined
FORWARD_BRIDGE  = "GGAGTGAGTACGGTGTGC"
R1_TAIL_BRIDGE  = "CCATCCAGCATCCAACTC"
R2_TAIL_BRIDGE  = "GCACACCGTACTCACTCC"
LEADER_LEN      = 27               # 4 random + 4 barcode + 1 + 18 bridge

# Allele motifs (must match config/targets.json)
MOTIFS = {
    ("chr01", "wt"):      "TTGTTGCTATAACTGGTCAAGTGCCACGT",
    ("chr01", "desired"): "TTGTTGCTATAACTGGTCAAGTGTCACGT",
    ("chr11", "wt"):      "TTGTTGCTATAACCGGTCAAGTGCCACGT",
    ("chr11", "desired"): "TTGTTGCTATAACCGGTCAAGTGTCACGT",
}

PLAN = [
    ("chr01", "wt",      60),
    ("chr01", "desired", 40),
    ("chr11", "wt",      30),
    ("chr11", "desired", 70),
]
N_UNDETERMINED = 25
FLANK_LEN      = 45                # keeps post-trim reads above the 100 bp QC threshold

_COMP = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(seq):
    return seq.translate(_COMP)[::-1]


def build_pair(rng, motif, forward_barcode, with_tail_bridge):
    """Return (r1_seq, r2_seq) for one synthetic amplicon molecule."""
    def rand_seq(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    amplicon = rand_seq(FLANK_LEN) + motif + rand_seq(FLANK_LEN)
    r1 = rand_seq(4) + forward_barcode + "T" + FORWARD_BRIDGE + amplicon
    r2 = rand_seq(4) + REVERSE_BARCODE + "T" + FORWARD_BRIDGE + reverse_complement(amplicon)
    assert len(r1) - len(amplicon) == LEADER_LEN

    if with_tail_bridge:
        # Short molecule: the read runs through into the opposite bridge.
        r1 += R1_TAIL_BRIDGE + rand_seq(8)
        r2 += R2_TAIL_BRIDGE + rand_seq(8)
    return r1, r2


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)

    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)

    records = []
    n = 0
    for allele, kind, count in PLAN:
        for i in range(count):
            n += 1
            r1, r2 = build_pair(rng, MOTIFS[(allele, kind)], FORWARD_BARCODE,
                                with_tail_bridge=(i % 5 == 0))
            records.append((f"synth{n}_{allele}_{kind}", r1, r2))

    for _ in range(N_UNDETERMINED):
        n += 1
        r1, r2 = build_pair(rng, MOTIFS[("chr01", "wt")], UNKNOWN_BARCODE,
                            with_tail_bridge=False)
        records.append((f"synth{n}_unknown_barcode", r1, r2))

    rng.shuffle(records)

    r1_path = out_dir / "SYNTH_S1_L001_R1_001.fastq.gz"
    r2_path = out_dir / "SYNTH_S1_L001_R2_001.fastq.gz"
    with gzip.open(r1_path, "wt", encoding="utf-8") as o1, \
         gzip.open(r2_path, "wt", encoding="utf-8") as o2:
        for name, r1, r2 in records:
            o1.write(f"@{name} 1:N:0:1\n{r1}\n+\n{'I' * len(r1)}\n")
            o2.write(f"@{name} 2:N:0:1\n{r2}\n+\n{'I' * len(r2)}\n")

    assigned = sum(c for _, _, c in PLAN)
    print(f"[synthetic] {len(records)} read pairs "
          f"({assigned} assigned, {N_UNDETERMINED} undetermined)")
    print(f"[synthetic] {r1_path}")
    print(f"[synthetic] {r2_path}")


if __name__ == "__main__":
    main()
