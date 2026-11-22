"""Unit tests for the SNIGRON analysis scripts."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def load(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[path.stem] = module
    spec.loader.exec_module(module)
    return module


snig = load(ROOT / "bin" / "snigron_variants.py")


# ------------------------------------------------------------------ SNIGRON
def _write_snps(path: Path, rows):
    pd.DataFrame(rows).to_csv(path, sep="\t", header=False, index=False)


def test_snigron_acquired_and_lost(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    common = [100, "A", "G", 100, 5, 5, 1000, 1000, 1, 1, "chr", "q"]
    _write_snps(tmp_path / "ctrl.snps", [common, [700, "C", "T", 700, 5, 5, 1000, 1000, 1, 1, "chr", "q"]])
    _write_snps(tmp_path / "trt.snps", [common, [550, "G", ".", 550, 5, 5, 1000, 999, 1, 1, "chr", "q"]])
    (tmp_path / "r.bed").write_text("chr\t500\t600\tgeneY\t0\t+\n")
    (tmp_path / "c.csv").write_text("comparison,control,treated\np,ctrl,trt\n")

    snig.main(["--snps", "ctrl.snps", "trt.snps", "--regions", "r.bed",
               "--comparisons", "c.csv", "--flank", "50"])
    diffs = pd.read_csv("snigron_comparisons.tsv", sep="\t", keep_default_na=False)
    got = set(zip(diffs["change"], diffs["ref_pos"], diffs["type"], diffs["region"]))
    assert got == {("acquired", 550, "deletion", "geneY"), ("lost", 700, "snp", "")}


def test_snigron_empty_snps_file(tmp_path):
    empty = tmp_path / "s.snps"
    empty.write_text("")
    assert snig.read_snps(empty, "s").empty


def test_collapse_indels_merges_runs():
    df = pd.DataFrame({
        "sample": "s", "ref_contig": "c",
        "ref_pos": [10, 20, 21, 22, 40, 40],
        "ref_base": ["A", "G", "T", "C", ".", "."],
        "alt_base": ["T", ".", ".", ".", "A", "G"],
        "type": ["snp", "deletion", "deletion", "deletion", "insertion", "insertion"],
    })
    out = snig.collapse_indels(df)
    assert out[["ref_pos", "type", "length", "ref_base", "alt_base"]].values.tolist() == [
        [10, "snp", 1, "A", "T"], [20, "deletion", 3, "GTC", "."], [40, "insertion", 2, ".", "AG"]]


def test_snigron_rotated_insertions_are_shared(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _write_snps(tmp_path / "ctrl.snps", [[300, ".", "A", 300, 5, 5, 1000, 1000, 1, 1, "chr", "q"],
                                         [300, ".", "T", 300, 5, 5, 1000, 1000, 1, 1, "chr", "q"]])
    _write_snps(tmp_path / "trt.snps", [[300, ".", "T", 300, 5, 5, 1000, 1000, 1, 1, "chr", "q"],
                                        [300, ".", "A", 300, 5, 5, 1000, 1000, 1, 1, "chr", "q"]])
    (tmp_path / "r.bed").write_text("chr\t0\t1000\tg\n")
    (tmp_path / "c.csv").write_text("comparison,control,treated\np,ctrl,trt\n")
    snig.main(["--snps", "ctrl.snps", "trt.snps", "--regions", "r.bed", "--comparisons", "c.csv"])
    assert pd.read_csv("snigron_comparisons.tsv", sep="\t").empty
