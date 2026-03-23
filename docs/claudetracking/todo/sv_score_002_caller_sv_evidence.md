# SV-SCORE-002: Use mean_variant_specific_molecules for SV evidence in caller

**Epic**: SV-SCORE (docs/claudetracking/overallplans/SV-SCORE.md)
**Priority**: high
**Depends on**: SV-SCORE-001
**Status**: todo

## Goal

In `kam-call/src/caller.rs`, use `mean_variant_specific_molecules` as the
molecule count `k` for SV-type variants. Done looks like: `call_variant` uses
the new field for `LargeDeletion`, `TandemDuplication`, and `Inversion`; all
SNV/indel callers continue to use `min_molecules`; tests cover both paths.

## Steps

1. Read `kam-call/src/caller.rs` in full, focusing on `call_variant` and
   `assign_filter`.

2. In `call_variant`, replace:
   ```rust
   let k = alt_evidence.min_molecules;
   let total = ref_evidence.min_molecules + alt_evidence.min_molecules;
   ```
   with:
   ```rust
   let is_sv = matches!(
       variant_type,
       VariantType::LargeDeletion | VariantType::TandemDuplication | VariantType::Inversion
   );
   let k = if is_sv {
       alt_evidence.mean_variant_specific_molecules.round() as u32
   } else {
       alt_evidence.min_molecules
   };
   let total = if is_sv {
       let ref_k = ref_evidence.mean_molecules.round() as u32;
       ref_k + k
   } else {
       ref_evidence.min_molecules + alt_evidence.min_molecules
   };
   ```

3. Add or update tests:
   - A test where a large (≥50bp) inversion path has high
     `mean_variant_specific_molecules` and low `min_molecules`, and confirm
     the call uses the mean (NALT reflects the mean, not the min).
   - A test confirming SNVs are unchanged (still use `min_molecules`).

4. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- Using `mean_molecules` for the reference total is a reasonable approximation.
  The reference path molecules span the entire target window and do not have
  the min bottleneck.
- FP risk from error k-mers: error k-mers typically have `min_molecules=1` and
  `mean_variant_specific_molecules=0.0` (they are reference k-mers or
  near-reference). This change should not inflate FPs.
- `strand_bias_p` still uses `min_simplex_fwd`/`min_simplex_rev`. For SV paths
  this produces misleading p-values. Fixing that is SV-STRAND epic.
