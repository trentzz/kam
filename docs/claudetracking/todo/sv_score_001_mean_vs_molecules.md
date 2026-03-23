# SV-SCORE-001: Add mean_variant_specific_molecules to PathEvidence

**Epic**: SV-SCORE (docs/claudetracking/overallplans/SV-SCORE.md)
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Add a `mean_variant_specific_molecules: f32` field to `PathEvidence` in
`kam-pathfind/src/score.rs` and compute it in `score_and_rank_paths`. Done
looks like: the field is present and populated for all alt paths; the reference
path always has value 0.0; existing tests continue to pass.

## Steps

1. Read `kam-pathfind/src/score.rs` in full.

2. Add the field to `PathEvidence`:
   ```rust
   /// Mean `n_molecules` across variant-specific k-mers only.
   ///
   /// Variant-specific k-mers are those not shared with the reference path.
   /// For large SV paths, this is a better evidence estimate than
   /// `min_molecules` because the minimum over 70–150 k-mer positions
   /// bottlenecks at 1–3 even at moderate VAF.
   ///
   /// Set to 0.0 for the reference path. Computed by `score_and_rank_paths`.
   pub mean_variant_specific_molecules: f32,
   ```

3. In `score_and_rank_paths`, after the `ref_kmer_set` is built:
   - For each alt path, collect `n_molecules` for all k-mers NOT in
     `ref_kmer_set`.
   - If the variant-specific set is empty, set `mean_variant_specific_molecules
     = alt_evidence.mean_molecules` (fall back to overall mean).
   - Otherwise compute the mean and set the field.
   - For the reference path, set the field to 0.0.

4. Update any struct initialisers or tests that construct `PathEvidence`
   directly (add `mean_variant_specific_molecules: 0.0` default).

5. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- The existing `min_variant_specific_duplex` field (same approach, different
  metric) is the model to follow. Use the same `ref_kmer_set` logic.
- Do not change how `min_molecules` is computed — it is still used by SNV/indel
  callers.
- No CLI changes in this task; those are in SV-SCORE-002.
