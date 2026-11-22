#!/usr/bin/env python3
"""Annotate MUMmer/dnadiff variants and compare samples against their controls.

Every sample is aligned to the same reference, so all variants share reference
coordinates and can be annotated with one regions BED and compared directly:
a variant is "acquired" in a treated sample when it is absent from its control
(same position and alleles), and "lost" when only the control carries it.

Inputs are the `.snps` files written by `dnadiff` (show-snps -rlTHC format).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SNPS_COLUMNS = ["ref_pos", "ref_base", "alt_base", "query_pos", "buff", "dist",
                "ref_len", "query_len", "ref_frame", "query_frame",
                "ref_contig", "query_contig"]
# SNPs are matched on their alleles. Indels are matched on position, type and length:
# the same insertion can be written with different (rotated) bases by MUMmer when the
# flanking sequence is repetitive, which would otherwise create false acquired/lost pairs.
KEY = ["ref_contig", "ref_pos", "type", "length", "match_allele"]


def read_snps(path: Path, sample: str) -> pd.DataFrame:
    """Read a dnadiff .snps file; an empty file yields an empty table."""
    if path.stat().st_size == 0:
        return pd.DataFrame(columns=SNPS_COLUMNS + ["sample", "type"])
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    df = df.iloc[:, :len(SNPS_COLUMNS)]
    df.columns = SNPS_COLUMNS[: df.shape[1]]
    df["ref_pos"] = df["ref_pos"].astype(int)
    df["sample"] = sample
    df["type"] = np.select(
        [df["ref_base"] == ".", df["alt_base"] == "."],
        ["insertion", "deletion"], default="snp")
    return df


def collapse_indels(df: pd.DataFrame) -> pd.DataFrame:
    """Merge the one-base rows dnadiff writes for each indel into single events.

    Deletions are runs of consecutive reference positions; insertions are runs of
    rows at the same reference position. `length` is the indel size (1 for SNPs).
    """
    if df.empty:
        return df.assign(length=pd.Series(dtype=int))
    df = df.sort_values(["sample", "ref_contig", "ref_pos"]).reset_index(drop=True)
    prev = df.shift()
    same_event = (
        (df["type"] != "snp")
        & (df["type"] == prev["type"])
        & (df["sample"] == prev["sample"])
        & (df["ref_contig"] == prev["ref_contig"])
        & np.where(df["type"] == "deletion",
                   df["ref_pos"] == prev["ref_pos"] + 1,
                   df["ref_pos"] == prev["ref_pos"])
    )
    df["event"] = (~same_event).cumsum()
    agg = df.groupby("event").agg(
        sample=("sample", "first"), ref_contig=("ref_contig", "first"),
        ref_pos=("ref_pos", "first"), type=("type", "first"),
        ref_base=("ref_base", "".join), alt_base=("alt_base", "".join),
        length=("ref_pos", "size"))
    agg["ref_base"] = agg["ref_base"].str.replace(".", "", regex=False).replace("", ".")
    agg["alt_base"] = agg["alt_base"].str.replace(".", "", regex=False).replace("", ".")
    return agg.reset_index(drop=True)


def read_regions(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str)
    regions = pd.DataFrame({"chrom": raw[0], "start": raw[1].astype(int), "end": raw[2].astype(int)})
    regions["region"] = raw[3] if raw.shape[1] > 3 else regions["chrom"] + ":" + raw[1] + "-" + raw[2]
    regions["region_strand"] = raw[5] if raw.shape[1] > 5 else "."
    return regions


def annotate(variants: pd.DataFrame, regions: pd.DataFrame, flank: int) -> pd.DataFrame:
    """Add a `region` column ('' when outside all regions; ';'-joined when several)."""
    labels = [[] for _ in range(len(variants))]
    pos = variants["ref_pos"].to_numpy()  # 1-based
    contig = variants["ref_contig"].to_numpy()
    for reg in regions.itertuples(index=False):
        if reg.region_strand == "+":
            lo, hi = reg.start - flank, reg.end
        elif reg.region_strand == "-":
            lo, hi = reg.start, reg.end + flank
        else:
            lo, hi = reg.start - flank, reg.end + flank
        # BED is 0-based half-open; dnadiff positions are 1-based
        hit = (contig == reg.chrom) & (pos > lo) & (pos <= hi)
        body = (pos > reg.start) & (pos <= reg.end)
        for i in np.flatnonzero(hit):
            labels[i].append(reg.region if body[i] else f"{reg.region}(flank)")
    out = variants.copy()
    out["region"] = [";".join(x) for x in labels]
    return out


def compare(variants: pd.DataFrame, control: str, treated: str, name: str) -> pd.DataFrame:
    variants = variants.assign(match_allele=np.where(
        variants["type"] == "snp", variants["ref_base"] + ">" + variants["alt_base"], ""))
    ctrl = variants[variants["sample"] == control]
    trt = variants[variants["sample"] == treated]
    acquired = trt.merge(ctrl[KEY], on=KEY, how="left", indicator=True)
    acquired = acquired[acquired["_merge"] == "left_only"].assign(change="acquired")
    lost = ctrl.merge(trt[KEY], on=KEY, how="left", indicator=True)
    lost = lost[lost["_merge"] == "left_only"].assign(change="lost")
    out = pd.concat([acquired, lost], ignore_index=True).drop(columns=["_merge", "match_allele"])
    out.insert(0, "comparison", name)
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--snps", nargs="+", type=Path, required=True,
                   help="dnadiff .snps files named <sample>.snps")
    p.add_argument("--regions", type=Path, required=True)
    p.add_argument("--comparisons", type=Path, required=True,
                   help="CSV with columns comparison,control,treated")
    p.add_argument("--flank", type=int, default=200)
    args = p.parse_args(argv)

    variants = pd.concat([read_snps(f, f.name.removesuffix(".snps")) for f in args.snps],
                         ignore_index=True)
    variants = collapse_indels(variants)
    variants = annotate(variants, read_regions(args.regions), args.flank)
    cols = ["sample", "ref_contig", "ref_pos", "ref_base", "alt_base", "type", "length", "region"]
    variants[cols].to_csv("snigron_variants_vs_reference.tsv", sep="\t", index=False)

    comps = pd.read_csv(args.comparisons)
    diffs = [compare(variants, r.control.strip(), r.treated.strip(), r.comparison.strip())
             for r in comps.itertuples(index=False)]
    diff_cols = ["comparison", "change"] + cols
    diffs = pd.concat(diffs, ignore_index=True)[diff_cols] if diffs else pd.DataFrame(columns=diff_cols)
    diffs.to_csv("snigron_comparisons.tsv", sep="\t", index=False)

    per_sample = (variants.groupby(["sample", "type"]).size().unstack(fill_value=0)
                  .reindex(columns=["snp", "insertion", "deletion"], fill_value=0))
    per_sample["in_regions"] = variants[variants["region"] != ""].groupby("sample").size()
    per_sample = per_sample.fillna(0).astype(int).reset_index()
    per_comp = (diffs.groupby(["comparison", "change", "type"]).size()
                .rename("n").reset_index())
    in_regions = diffs[diffs["region"] != ""]

    per_sample.to_csv("snigron_summary.tsv", sep="\t", index=False)
    with pd.ExcelWriter("snigron_report.xlsx") as xl:
        per_sample.to_excel(xl, sheet_name="per_sample", index=False)
        per_comp.to_excel(xl, sheet_name="per_comparison", index=False)
        in_regions.to_excel(xl, sheet_name="changes_in_regions", index=False)
        variants[variants["region"] != ""][cols].to_excel(
            xl, sheet_name="variants_in_regions", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
