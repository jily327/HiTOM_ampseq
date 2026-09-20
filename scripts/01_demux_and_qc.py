#!/usr/bin/env python3
"""
01_demux_and_qc.py  --  Barcode demultiplexing and quality control.

Faithful port of the validated pipeline:
  - demux_hitom.py   : barcode extraction from R1/R2 + 27-bp leader trim
  - trim_tail_bridges.py : clip tail-bridge sequences from trimmed reads

Additional QC (new, not in validated motif-count path):
  - MD5 check of input FASTQs against .md5 sidecar files if present
  - Read-count and read-length summary
  - Barcode assignment rate summary
  - Resolved barcode map written to output

Outputs (all under <output>/01_demux/):
  demux_trimmed/         per-sample FASTQs after barcode demux + 27-bp trim
  demux_clean/           per-sample FASTQs after additional tail-bridge clip
  demux_counts.tsv       read pairs per sample (including undetermined)
  tail_trim_summary.tsv  tail-bridge clip counts per sample
  barcode_assignment_summary.csv  overall assignment rate
  read_length_summary.csv         post-trim read length stats per sample
  barcode_map_resolved.csv        barcode names → sequences used this run
  md5_summary.txt                 MD5 check results
"""

import argparse
import gzip
import hashlib
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

sys.path.insert(0, str(Path(__file__).resolve().parent))
import utils

PIPELINE_ROOT    = Path(__file__).resolve().parent.parent
DEFAULT_METADATA = PIPELINE_ROOT / "config" / "sample_metadata.csv"
DEFAULT_PROTOCOL = PIPELINE_ROOT / "config" / "hitom_protocol_constants.json"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def md5_file(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input",    required=True,
                   help="Directory containing raw paired FASTQ.gz files")
    p.add_argument("--output",   required=True,
                   help="Root output directory (outputs go into <output>/01_demux/)")
    p.add_argument("--metadata", default=str(DEFAULT_METADATA),
                   help=f"Sample metadata CSV  [default: pipeline/config/sample_metadata.csv]")
    p.add_argument("--protocol-constants", dest="protocol_constants",
                   default=str(DEFAULT_PROTOCOL),
                   help=f"Hi-TOM protocol constants JSON  [default: pipeline/config/hitom_protocol_constants.json]")
    # Optional overrides for protocol parameters
    p.add_argument("--r1-bc-start",   type=int, default=None, dest="r1_bc_start")
    p.add_argument("--r1-bc-end",     type=int, default=None, dest="r1_bc_end")
    p.add_argument("--r2-bc-start",   type=int, default=None, dest="r2_bc_start")
    p.add_argument("--r2-bc-end",     type=int, default=None, dest="r2_bc_end")
    p.add_argument("--trim-start-r1", type=int, default=None, dest="trim_start_r1")
    p.add_argument("--trim-start-r2", type=int, default=None, dest="trim_start_r2")
    p.add_argument("--r1-tail-bridge", default=None, dest="r1_tail_bridge")
    p.add_argument("--r2-tail-bridge", default=None, dest="r2_tail_bridge")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # --- Load protocol constants ---
    proto   = utils.load_json(args.protocol_constants)
    fw_dict = {k.upper(): v.upper() for k, v in proto["forward_barcodes"].items()
               if not k.startswith("_")}
    rv_dict = {k.upper(): v.upper() for k, v in proto["reverse_barcodes"].items()
               if not k.startswith("_")}
    defs    = proto["defaults"]

    # Parameters: CLI override → protocol constant default
    r1_bc_start    = args.r1_bc_start    if args.r1_bc_start    is not None else defs["r1_bc_start"]
    r1_bc_end      = args.r1_bc_end      if args.r1_bc_end      is not None else defs["r1_bc_end"]
    r2_bc_start    = args.r2_bc_start    if args.r2_bc_start    is not None else defs["r2_bc_start"]
    r2_bc_end      = args.r2_bc_end      if args.r2_bc_end      is not None else defs["r2_bc_end"]
    trim_r1        = args.trim_start_r1  if args.trim_start_r1  is not None else defs["trim_start_r1"]
    trim_r2        = args.trim_start_r2  if args.trim_start_r2  is not None else defs["trim_start_r2"]
    r1_tail_bridge = (args.r1_tail_bridge or proto["tail_bridges"]["r1_tail_bridge"]).upper()
    r2_tail_bridge = (args.r2_tail_bridge or proto["tail_bridges"]["r2_tail_bridge"]).upper()

    print(f"[01_demux] Protocol parameters:")
    print(f"  Barcode positions  R1[{r1_bc_start}:{r1_bc_end}]  R2[{r2_bc_start}:{r2_bc_end}]")
    print(f"  Leader trim        R1={trim_r1} bp  R2={trim_r2} bp")
    print(f"  Tail bridges       R1={r1_tail_bridge}  R2={r2_tail_bridge}")

    # --- Load metadata + build barcode lookup ---
    meta_rows, bc_lookup = utils.load_metadata(args.metadata, fw_dict, rv_dict)
    sample_ids = [r["sample_id"] for r in meta_rows]
    print(f"[01_demux] Loaded {len(sample_ids)} samples from metadata")

    # --- Locate input FASTQs ---
    in_dir = Path(args.input)
    r1_path, r2_path = utils.locate_fastq_pair(in_dir)
    print(f"[01_demux] R1: {r1_path}")
    print(f"[01_demux] R2: {r2_path}")

    # --- Output directories ---
    out_dir     = Path(args.output) / "01_demux"
    trimmed_dir = out_dir / "demux_trimmed"
    clean_dir   = out_dir / "demux_clean"
    out_dir.mkdir(parents=True, exist_ok=True)
    trimmed_dir.mkdir(exist_ok=True)
    clean_dir.mkdir(exist_ok=True)

    # --- MD5 check ---
    print("[01_demux] Computing MD5 checksums...")
    r1_md5 = md5_file(r1_path)
    r2_md5 = md5_file(r2_path)
    md5_lines = [
        f"R1: {r1_path.name}  md5={r1_md5}",
        f"R2: {r2_path.name}  md5={r2_md5}",
    ]
    for fpath, fmd5 in [(r1_path, r1_md5), (r2_path, r2_md5)]:
        sidecar = Path(str(fpath) + ".md5")
        if sidecar.exists():
            expected = sidecar.read_text().split()[0].strip()
            status = "MATCH" if fmd5 == expected else "MISMATCH"
            md5_lines.append(f"  {fpath.name}: {status}  expected={expected}")
            print(f"[01_demux]   {fpath.name}: {status}")
    (out_dir / "md5_summary.txt").write_text("\n".join(md5_lines) + "\n")

    # -----------------------------------------------------------------------
    # PASS 1 — Demultiplex + 27-bp leader trim
    # Faithful port of demux_hitom.py
    # -----------------------------------------------------------------------
    print("[01_demux] Pass 1: demultiplexing + leader trim...")
    handles_trimmed = {}
    counts   = defaultdict(int)
    lengths  = defaultdict(list)   # sample → [post-trim R1 lengths]

    def get_trimmed_handles(sample):
        if sample not in handles_trimmed:
            h1 = gzip.open(trimmed_dir / f"{sample}_R1.fastq.gz", "wt")
            h2 = gzip.open(trimmed_dir / f"{sample}_R2.fastq.gz", "wt")
            handles_trimmed[sample] = (h1, h2)
        return handles_trimmed[sample]

    with utils.open_gz(r1_path) as f1, utils.open_gz(r2_path) as f2:
        for (h1, s1, p1, q1), (h2, s2, p2, q2) in zip(
                utils.fastq_iter(f1), utils.fastq_iter(f2)):

            fb = s1.upper()[r1_bc_start:r1_bc_end]
            rb = s2.upper()[r2_bc_start:r2_bc_end]
            sample = bc_lookup.get((fb, rb), "undetermined")

            # Trim 27-bp leader
            s1t, q1t = s1[trim_r1:], q1[trim_r1:]
            s2t, q2t = s2[trim_r2:], q2[trim_r2:]

            o1, o2 = get_trimmed_handles(sample)
            o1.write(f"{h1}\n{s1t}\n{p1}\n{q1t}\n")
            o2.write(f"{h2}\n{s2t}\n{p2}\n{q2t}\n")
            counts[sample] += 1
            if sample != "undetermined":
                lengths[sample].append(len(s1t))

    for h1, h2 in handles_trimmed.values():
        h1.close(); h2.close()

    total_reads    = sum(counts.values())
    assigned_reads = sum(v for k, v in counts.items() if k != "undetermined")
    assign_rate    = assigned_reads / total_reads * 100 if total_reads else 0
    print(f"[01_demux]   Total read pairs : {total_reads:,}")
    print(f"[01_demux]   Assigned         : {assigned_reads:,}  ({assign_rate:.1f}%)")
    print(f"[01_demux]   Undetermined     : {counts['undetermined']:,}")

    # Write demux_counts.tsv
    utils.write_tsv(
        [{"sample": s, "read_pairs": c} for s, c in sorted(counts.items())],
        out_dir / "demux_counts.tsv",
    )

    # Write barcode_assignment_summary.csv
    utils.write_csv([{
        "total_read_pairs":   total_reads,
        "assigned_read_pairs": assigned_reads,
        "undetermined_read_pairs": counts["undetermined"],
        "assignment_rate_pct": f"{assign_rate:.2f}",
        "n_samples_with_reads": len([s for s in counts if s != "undetermined"]),
    }], out_dir / "barcode_assignment_summary.csv")

    # Write read_length_summary.csv
    len_rows = []
    for sid in sorted(lengths):
        ls = lengths[sid]
        if ls:
            len_rows.append({
                "sample":          sid,
                "n_reads":         len(ls),
                "mean_R1_len_bp":  f"{mean(ls):.1f}",
                "median_R1_len_bp": f"{median(ls):.0f}",
                "min_R1_len_bp":   min(ls),
                "max_R1_len_bp":   max(ls),
            })
    utils.write_csv(len_rows, out_dir / "read_length_summary.csv")

    # -----------------------------------------------------------------------
    # PASS 2 — Tail-bridge clipping (demux_trimmed → demux_clean)
    # Faithful port of trim_tail_bridges.py
    # -----------------------------------------------------------------------
    print("[01_demux] Pass 2: tail-bridge clipping...")
    tail_summary = []

    for r1_trim in sorted(trimmed_dir.glob("*_R1.fastq.gz")):
        sid     = r1_trim.name.replace("_R1.fastq.gz", "")
        r2_trim = trimmed_dir / f"{sid}_R2.fastq.gz"
        if not r2_trim.exists():
            continue

        out_r1 = clean_dir / r1_trim.name
        out_r2 = clean_dir / r2_trim.name

        total = r1_clipped = r2_clipped = 0

        with utils.open_gz(r1_trim) as f1, utils.open_gz(r2_trim) as f2, \
             gzip.open(out_r1, "wt") as o1, gzip.open(out_r2, "wt") as o2:

            for (h1, s1, p1, q1), (h2, s2, p2, q2) in zip(
                    utils.fastq_iter(f1), utils.fastq_iter(f2)):

                idx1 = s1.upper().find(r1_tail_bridge)
                if idx1 >= 0:
                    s1, q1 = s1[:idx1], q1[:idx1]
                    r1_clipped += 1

                idx2 = s2.upper().find(r2_tail_bridge)
                if idx2 >= 0:
                    s2, q2 = s2[:idx2], q2[:idx2]
                    r2_clipped += 1

                o1.write(f"{h1}\n{s1}\n{p1}\n{q1}\n")
                o2.write(f"{h2}\n{s2}\n{p2}\n{q2}\n")
                total += 1

        tail_summary.append({
            "sample":                  sid,
            "total_read_pairs":        total,
            "R1_tail_bridge_trimmed":  r1_clipped,
            "R2_tail_bridge_trimmed":  r2_clipped,
        })

    utils.write_tsv(tail_summary, out_dir / "tail_trim_summary.tsv")

    # -----------------------------------------------------------------------
    # Write resolved barcode map
    # -----------------------------------------------------------------------
    bc_map_rows = []
    for row in meta_rows:
        fw_seq = utils.resolve_barcode_seq(row["forward_barcode"], fw_dict, rv_dict, "forward")
        rv_seq = utils.resolve_barcode_seq(row["reverse_barcode"], fw_dict, rv_dict, "reverse")
        bc_map_rows.append({
            "sample_id":             row["sample_id"],
            "forward_barcode_name":  row["forward_barcode"],
            "forward_barcode_seq":   fw_seq,
            "reverse_barcode_name":  row["reverse_barcode"],
            "reverse_barcode_seq":   rv_seq,
        })
    utils.write_csv(bc_map_rows, out_dir / "barcode_map_resolved.csv")

    print(f"[01_demux] Done. All outputs in {out_dir}")


if __name__ == "__main__":
    main()
