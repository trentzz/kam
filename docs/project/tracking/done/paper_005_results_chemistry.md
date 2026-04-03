# PAPER-005: Results Section — Cross-Chemistry Generalisation

**Epic**: PAPER (docs/claudetracking/overallplans/PAPER.md)
**Priority**: high
**Depends on**: pub_005_cross_dataset_table.md
**Status**: todo

## Goal

Write the cross-chemistry results subsection: public dataset performance
demonstrating that kam works beyond Twist UMI duplex chemistry. Done looks
like: a complete draft (≈350 words) committed to
`docs/paper/sections/results_chemistry.md`.

## Steps

1. Read `docs/benchmarking/public/cross_dataset_table.tsv` and the per-dataset
   results notes from PUB-002, PUB-003, and PUB-004.

2. Write `docs/paper/sections/results_chemistry.md` covering:

   **UMI Clustering Benchmark (≈100 words)**
   - Chemistry: simplex 12 bp UMI, no skip.
   - Molecule assembly results: molecule count, UMI collision rate.
   - Comparison to published HUMID results for the same sample.
   - Key message: kam correctly handles non-Twist UMI chemistry.

   **SEQC2 HCC1395 (≈100 words)**
   - SNV/indel sensitivity and precision against Tier 1 truth.
   - Mode: tumour-informed only (or discovery if no truth set available).
   - Caveat: dataset may not have UMIs; note if a bypass mode was used.

   **COLO829 SV (≈100 words)**
   - SV sensitivity across the 68 validated structural variants.
   - Sensitivity per SV type.
   - Caveat: WGS, no UMI; bypass mode used.

3. If any public benchmark is incomplete at the time of writing, note it
   with `[TODO: run pending]` rather than omitting the section.

## Notes

- The cross-chemistry results support the claim that kam is not Twist-specific.
  Frame each dataset as testing a different aspect of generalisation.
- Write limitations honestly: "Because SRR6794144 does not include a variant
  truth set, we evaluated molecule assembly quality rather than variant
  detection accuracy."
