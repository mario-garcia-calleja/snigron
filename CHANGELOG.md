# Changelog

## 2.0.0

Complete rewrite of the Bachelor's thesis version (v1).

- **Comparisons in reference coordinates.** v1 aligned treated vs control directly and annotated
  those positions with reference gene coordinates, which is only valid if the control assembly is
  collinear with the reference. v2 aligns every sample to the reference and derives acquired/lost
  variants by comparing call sets.
- **Indels as events.** `dnadiff` writes one row per indel base; v2 collapses them into single
  events and matches indels by position and length (MUMmer can write the same insertion with
  rotated bases in repetitive context).
- **Promoters.** Variants upstream of a region (`--flank`, strand-aware) are reported.
- IntegronFinder functional annotation is optional (`--func_hmm_dir`); `--skip_integrons` added.
- Samplesheet-driven, organism-agnostic, resource limits and execution profiles.
- Analysis script rewritten with argparse; TSV outputs plus one Excel report; unit tests;
  synthetic test dataset and GitHub Actions CI.
