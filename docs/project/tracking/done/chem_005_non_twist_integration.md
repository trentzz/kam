# CHEM-005: Integration Test with Non-Twist Chemistry

**Epic**: CHEM-CONFIG (docs/claudetracking/overallplans/CHEM-CONFIG.md)
**Priority**: critical
**Depends on**: chem_003_generalise_clustering.md, chem_004_read_structure_presets.md
**Status**: todo

## Goal

Write an end-to-end integration test that runs the full kam pipeline on
synthetic data with a 12 bp UMI and no skip. Done looks like: the test
generates synthetic FASTQ with 12 bp UMI structure, runs `kam run` with
`preset = simplex-12bp`, and produces a calls VCF with the expected SNV
detected.

## Steps

1. Write a FASTQ generator function (or script) that produces paired-end reads
   with a 12 bp random UMI prefix and a 2 × 100 bp read body containing one
   SNV at a known position. Use a fixed seed for determinism.

2. Write a test config TOML:
   ```toml
   [chemistry]
   preset = "simplex-12bp"

   [calling]
   min_confidence = 0.95
   min_alt_molecules = 1
   ```

3. Write the integration test in `kam/tests/` (or a new file
   `kam/tests/chem_integration.rs`):
   - Generate the synthetic FASTQs into a `tempdir`.
   - Run `kam run --config test.toml --r1 ... --r2 ... --targets ...`.
   - Assert the output VCF contains the expected SNV with PASS filter.
   - Assert no other variants are called (no false positives from UMI parsing
     with wrong length).

4. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

5. Document the test design in a comment at the top of the test file: what
   chemistry is being simulated, what the expected output is, and why this
   test is important.

## Notes

- The synthetic FASTQ does not need a skip field (skip_length = 0). The
  first 12 bytes of each read are the UMI; byte 13 onwards is the template.
- If varforge does not support 12 bp UMI, write the FASTQ generator in Rust
  directly in the test file.
- This test guards against regressions where Twist-specific assumptions
  silently break non-Twist input.
