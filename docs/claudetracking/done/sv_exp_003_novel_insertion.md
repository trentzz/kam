# SV-EXP-003: NovelInsertion Classification

**Epic**: SV-EXPAND (docs/claudetracking/overallplans/SV-EXPAND.md)
**Priority**: high
**Depends on**: sv_exp_001_new_variant_types.md
**Status**: todo

## Goal

Classify a path as `NovelInsertion` instead of `TandemDuplication` when the
inserted sequence (≥50 bp) has no tandem repeat relationship to the nearby
reference. Done looks like: a synthetic insertion of non-repetitive sequence is
classified as `NovelInsertion`; a tandem duplication of a nearby reference
window is still classified as `TandemDuplication`.

## Steps

1. Read `kam-call/src/classify.rs` in full, focusing on the
   `TandemDuplication` classification branch.

2. Add a helper function `is_tandem_repeat(inserted: &[u8], ref_window: &[u8])
   -> bool`:
   - Check whether the inserted sequence is a copy (or near-copy, ≤2 mismatches)
     of any substring of `ref_window` of the same length.
   - Return `true` if a tandem repeat relationship is found.

3. In the `TandemDuplication` branch, after the insertion is detected (alt
   longer than ref by ≥50 bp):
   - Extract the inserted bases from the alt path sequence.
   - Extract a reference window of 2× the insertion length centred on the
     insertion site.
   - Call `is_tandem_repeat(inserted, ref_window)`.
   - If `false`, classify as `NovelInsertion` instead.

4. Write unit tests:
   - Insertion that duplicates a nearby ref segment → `TandemDuplication`.
   - Insertion of random non-reference sequence (≥50 bp) → `NovelInsertion`.
   - Insertion of 49 bp → neither (below threshold), stays as existing behaviour.
   - Near-copy with 1 mismatch → `TandemDuplication` (within tolerance).
   - Near-copy with 3 mismatches → `NovelInsertion` (exceeds tolerance).

5. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- The reference window for the tandem repeat check must be accessible from the
  classifier. If only path sequences (not raw reference bytes) are available,
  use the ref path sequence as the window.
- Homopolymer runs (e.g. AAAAAAA...A, 50 bp) should remain `TandemDuplication`.
  The `is_tandem_repeat` check will naturally handle this because the repeat
  unit is a substring of the reference.
- The 50 bp insertion threshold and 2-mismatch tolerance are named constants.
