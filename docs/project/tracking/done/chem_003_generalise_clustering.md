# CHEM-003: Generalise Clustering and Consensus for Arbitrary UMI Lengths

**Epic**: CHEM-CONFIG (docs/claudetracking/overallplans/CHEM-CONFIG.md)
**Priority**: critical
**Depends on**: chem_002_generalise_parser.md
**Status**: todo

## Goal

Ensure Hamming clustering and consensus calling work correctly with any UMI
length. Remove any UMI-length assertions or compile-time assumptions in
`kam-assemble/src/clustering.rs` and `kam-assemble/src/consensus.rs`. Done
looks like: clustering and consensus produce correct results for both 5 bp and
12 bp UMIs in unit tests, with no hardcoded length constants remaining.

## Steps

1. Read `kam-assemble/src/clustering.rs` and `kam-assemble/src/consensus.rs`
   in full.

2. Search for any assertions or constants that assume a specific UMI length:
   - `assert_eq!(umi.len(), 5)` or similar
   - Array indexing `umi[..5]`
   - Hardcoded loop bounds based on UMI length

3. The Hamming distance function likely operates on byte slices already. Verify
   it uses `a.len()` (or checks `a.len() == b.len()`) rather than a constant.
   If it asserts equal length, keep that assertion (it is correct behaviour).

4. In consensus calling, check whether the UMI is referenced. If consensus
   only operates on read bases (not the UMI itself), no change is needed.
   If it uses UMI length for indexing into reads, update those offsets to use
   the configured `umi_length + skip_length`.

5. Add or update tests:
   - Hamming clustering with 12 bp UMIs: two UMIs differing at position 11
     are correctly grouped.
   - Clustering correctly rejects UMIs with more than `max_hamming` differences
     at any length.

6. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- If `clustering.rs` already uses slice operations throughout and has no
  length assumptions, the main work is just verification and adding the 12 bp
  test.
- `consensus.rs` may reference read offsets like `read[umi_length + skip_length..]`
  to skip the UMI/skip prefix when building consensus. These must use the
  configured lengths, not constants.
