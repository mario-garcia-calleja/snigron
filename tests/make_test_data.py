#!/usr/bin/env python3
"""Generate the small synthetic dataset used by `-profile test` and the unit tests.

Everything is simulated (no real sequencing data is distributed with the repository):
  reference.fasta  20 kb random genome ('chrT')
  regions.bed      three 'genes' on chrT
  *.fasta          assemblies derived from the reference with known SNPs and a deletion

Run from the repository root:  python tests/make_test_data.py
"""
from __future__ import annotations

import random
from pathlib import Path

OUT = Path(__file__).parent / "data"
CHROM = "chrT"
SEED = 42


def make_reference(rng: random.Random, length: int = 20_000) -> str:
    seq = list(rng.choices("ACGT", k=length))
    for pos in range(100, length - 100, 250):  # plant a GATC every 250 bp
        seq[pos:pos + 4] = list("GATC")
    return "".join(seq)


def write_fasta(path: Path, name: str, seq: str) -> None:
    lines = [seq[i:i + 70] for i in range(0, len(seq), 70)]
    path.write_text(f">{name}\n" + "\n".join(lines) + "\n")


def mutate(seq: str, snps: dict[int, str], deletions: list[tuple[int, int]]) -> str:
    s = list(seq)
    for pos, alt in snps.items():
        s[pos] = alt
    for start, length in sorted(deletions, reverse=True):
        del s[start:start + length]
    return "".join(s)


def main() -> None:
    rng = random.Random(SEED)
    OUT.mkdir(parents=True, exist_ok=True)

    ref = make_reference(rng)
    write_fasta(OUT / "reference.fasta", CHROM, ref)
    (OUT / "regions.bed").write_text(
        f"{CHROM}\t2000\t4000\tgene_A\t0\t+\n"
        f"{CHROM}\t9000\t11000\tgene_B\t0\t-\n"
        f"{CHROM}\t15000\t16500\tgene_C\t0\t.\n")

    # SNIGRON: control carries one background SNP, treated adds SNPs in genes + a deletion
    background = {5000: "T" if ref[5000] != "T" else "G"}
    acquired = {3000: "A" if ref[3000] != "A" else "C",      # gene_A
                10_000: "G" if ref[10_000] != "G" else "T"}  # gene_B
    write_fasta(OUT / "ctrl_1.fasta", "ctrl_1", mutate(ref, background, []))
    write_fasta(OUT / "trt_1.fasta", "trt_1",
                mutate(ref, {**background, **acquired}, [(15_500, 12)]))
    (OUT / "samplesheet.csv").write_text(
        "sample,fasta\nctrl_1,ctrl_1.fasta\ntrt_1,trt_1.fasta\n")
    (OUT / "comparisons.csv").write_text("comparison,control,treated\npair1,ctrl_1,trt_1\n")


if __name__ == "__main__":
    main()
