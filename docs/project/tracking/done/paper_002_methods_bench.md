# PAPER-002: Methods Section — Benchmarking Design

**Epic**: PAPER (docs/claudetracking/overallplans/PAPER.md)
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Write the benchmarking design subsection of the methods: datasets, evaluation
metrics, and comparison strategy. Done looks like: a complete draft (≈500
words) committed to `docs/paper/sections/methods_bench.md`.

## Steps

1. Read `docs/benchmarking/` directory structure and the BENCH-VARFORGE,
   BENCH-SV-NEW, PUB-BENCH, and ALIGN-COMPARE epic plans to understand the
   full benchmarking scope.

2. Write `docs/paper/sections/methods_bench.md` covering:

   **Synthetic benchmarks (≈150 words)**
   - varforge-generated datasets.
   - 25 VAF levels × 2 replicates per variant class.
   - Variant classes: SNV, indel, DEL, DUP, INV, InvDel, NovelInsertion, Fusion.
   - Both discovery and tumour-informed modes evaluated.

   **Public datasets (≈150 words)**
   - UMI Clustering Benchmark (SRR6794144): chemistry generalisation.
   - SEQC2 HCC1395: SNV/indel accuracy on real data.
   - COLO829: SV detection on 68 validated structural variants.
   - Data provenance and download accessions.

   **Alignment comparison (≈100 words)**
   - 24 titration samples, 375 variants.
   - Comparison method: per-variant concordance table.
   - Alignment pipeline: [RaSCALL / HUMID + km — confirm tool names from thesis].

   **Evaluation metrics (≈100 words)**
   - Sensitivity = TP / (TP + FN).
   - Precision = TP / (TP + FP) where a truth set allows it.
   - SV breakpoint matching tolerance: ±50 bp.
   - SNV/indel matching: exact chrom + pos + ref + alt.

3. Include `[TODO: fill accession numbers after download]` for dataset
   accessions not yet confirmed.

## Notes

- Write in Australian English, active voice, no em dashes.
- Methods sections are written in past tense ("we evaluated", "we ran").
- This section can be drafted before benchmarks are complete.
