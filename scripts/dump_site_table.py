#!/usr/bin/env python3
"""Project a finished MitoRiboPy run's BED + offset CSVs into a site
table that ``mitoribopy periodicity`` can re-score.

Reads:
    <run_root>/align/bed/<sample>.bed
    <run_root>/rpf/offset_diagnostics/csv/per_sample_offset/<sample>/offset_applied.csv

Writes:
    <output> (TSV) with columns:
        sample, gene, transcript_id, read_length, site_type,
        site_pos, cds_start, cds_end

The site_type column is fixed to ``p`` (P-site). To dump A-site rows,
re-run with ``--site a``; the helper just shifts the assigned
coordinate by +3 nt the same way the rest of the pipeline does.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from mitoribopy.analysis.periodicity import compute_p_site_positions
from mitoribopy.data.reference_data import load_annotation_table


def _read_bed(path: Path, sample: str) -> pd.DataFrame:
    """Read one BED6 file with the column names the pipeline uses."""
    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["chrom", "start", "end", "name", "score", "strand"],
        comment="#",
    )
    df["sample_name"] = sample
    df["read_length"] = (df["end"].astype(int) - df["start"].astype(int)).astype(int)
    return df


def _build_offset_map(per_sample_dir: Path, samples: list[str]) -> dict[str, pd.DataFrame]:
    out: dict[str, pd.DataFrame] = {}
    for s in samples:
        csv = per_sample_dir / s / "offset_applied.csv"
        if csv.is_file():
            out[s] = pd.read_csv(csv)
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--run-root", required=True, help="Path to the finished MitoRiboPy run root.")
    p.add_argument("--output",   required=True, help="Destination TSV path.")
    p.add_argument("--site", choices=("p", "a"), default="p")
    p.add_argument("--strain", default="h.sapiens")
    p.add_argument("--annotation-file", default=None,
                   help="Override built-in annotation (custom strains).")
    p.add_argument("--offset-type", choices=("5", "3"), default="5")
    args = p.parse_args()

    run_root = Path(args.run_root)
    bed_dir = run_root / "align" / "bed"
    offset_dir = run_root / "rpf" / "offset_diagnostics" / "csv" / "per_sample_offset"

    if not bed_dir.is_dir():
        print(f"ERROR: missing {bed_dir}", file=sys.stderr)
        return 2
    if not offset_dir.is_dir():
        print(f"ERROR: missing {offset_dir}", file=sys.stderr)
        return 2

    samples = sorted(p.stem for p in bed_dir.glob("*.bed"))
    if not samples:
        print(f"ERROR: no BED files under {bed_dir}", file=sys.stderr)
        return 2

    annotation_df = load_annotation_table(
        preset=args.strain, annotation_file=args.annotation_file,
    )
    chrom_to_tx = annotation_df.set_index("sequence_name")["transcript"].to_dict()
    tx_meta = annotation_df.set_index("transcript")[
        ["start_codon", "l_cds", "display_name"]
    ]

    offset_map = _build_offset_map(offset_dir, samples)
    if not offset_map:
        print(f"ERROR: no per-sample offset CSVs under {offset_dir}", file=sys.stderr)
        return 2

    rows: list[pd.DataFrame] = []
    for sample in samples:
        bed = _read_bed(bed_dir / f"{sample}.bed", sample)
        psite = compute_p_site_positions(
            bed,
            sample=sample,
            selected_offsets_by_sample=offset_map,
            selected_offsets_combined=None,
            offset_type=args.offset_type,
            offset_site="p",
        )
        if psite.empty:
            continue
        psite["transcript_id"] = (
            psite["chrom"].map(chrom_to_tx).fillna(psite["chrom"]).astype(str)
        )
        psite = psite[psite["transcript_id"].isin(tx_meta.index)]
        if psite.empty:
            continue
        psite["cds_start"] = psite["transcript_id"].map(tx_meta["start_codon"]).astype(int)
        psite["cds_end"] = (
            psite["cds_start"] + psite["transcript_id"].map(tx_meta["l_cds"]).astype(int)
        )
        psite["gene"] = psite["transcript_id"].map(tx_meta["display_name"]).fillna(
            psite["transcript_id"]
        )
        site_pos = psite["P_site"].astype(int)
        if args.site == "a":
            site_pos = site_pos + 3
        out = pd.DataFrame({
            "sample": sample,
            "gene": psite["gene"].astype(str),
            "transcript_id": psite["transcript_id"].astype(str),
            "read_length": psite["read_length"].astype(int),
            "site_type": args.site,
            "site_pos": site_pos.astype(int),
            "cds_start": psite["cds_start"].astype(int),
            "cds_end": psite["cds_end"].astype(int),
        })
        rows.append(out)

    if not rows:
        print("ERROR: no usable rows after P-site assignment", file=sys.stderr)
        return 2

    final = pd.concat(rows, ignore_index=True)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(out_path, sep="\t", index=False)
    print(
        f"[dump_site_table] wrote {len(final):,} rows across "
        f"{final['sample'].nunique()} sample(s) -> {out_path}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
