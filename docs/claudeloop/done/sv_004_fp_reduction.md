# Task SV-004: Add Monitoring Mode to SV Benchmark and Reduce FPs

## Problem

SV targets currently produce many SNV/MNV false positives from sequencing
error reads in the target window. At 5% VAF the DEL target alone has ~300
calls, of which only 1 is the true deletion. All others are low-molecule
error variants.

## What to do

1. Update `docs/benchmarking/sv/scripts/run_kam_sv.sh` to run a second pass
   per VAF level using `--target-variants truth_svs_vaf${TAG}.vcf`.
   Write results to `kam_vaf${TAG}_tumour_informed/variants.tsv`.

2. In tumour-informed mode, verify:
   - Only the true SV (DEL / DUP / INV) is labelled PASS.
   - All SNV/MNV false positives are relabelled NotTargeted.
   - FP count drops to 0.

3. Check whether the truth VCF symbolic allele format (`<DEL>`, `<DUP>`,
   `<INV>`) matches what kam outputs. The tumour-informed mode matching logic in
   `kam-call/src/targeting.rs` uses exact (CHROM, POS, REF, ALT) comparison.
   Symbolic alleles may require special handling.

   If symbolic allele matching does not work, investigate and fix the matching
   logic, or switch the truth VCF to sequence-resolved alleles for the SV
   benchmark.

## Definition of done

- `run_kam_sv.sh` runs both discovery and tumour-informed passes.
- Monitoring mode results in 0 FPs for each SV type at each VAF level.
- Discovery mode FP count is documented for comparison.
