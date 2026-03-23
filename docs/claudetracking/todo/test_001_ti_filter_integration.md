# TEST-001: Integration tests for tumour-informed filtering and multi-format output

**Epic**: none (standalone test work)
**Priority**: medium
**Depends on**: none
**Status**: todo

## Goal

Add integration tests covering:
1. The `--target-variants` tumour-informed filter path in both `call` and `run`
   subcommands.
2. The `--output-format tsv,vcf` multi-format output path.

Done looks like: both paths have at least one round-trip test that verifies
FP calls are filtered and the output files exist.

## Steps

1. Read `kam/src/commands/call.rs` tests and `kam/src/commands/run.rs` tests.

2. Add a test to `call.rs`:
   - Make a synthetic VCF with a known target variant.
   - Run `call_variants` with `--target-variants` pointing to the VCF.
   - Assert that a call at the target position PASSES and a call at a non-
     target position is `NotTargeted`.

3. Add a test to `run.rs`:
   - Use existing synthetic FASTQ infrastructure.
   - Run with `--target-variants` and assert only targeted calls appear.

4. Add a test for `--output-format tsv,vcf`:
   - Run with both formats.
   - Assert both output files exist and are non-empty.

5. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- From codebase review findings 6.1 and 6.3.
- Finding 6.2 (pathfind with variant detection) is a separate task since it
  requires a different test setup.
