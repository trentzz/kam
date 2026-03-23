# SV-STRAND-001: Add sv_strand_bias_threshold to CallerConfig

**Epic**: SV-STRAND (docs/claudetracking/overallplans/SV-STRAND.md)
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Add `sv_strand_bias_threshold: f64` to `CallerConfig` and use it in
`assign_filter` for SV-type variants. Done looks like: INV calls at ≥2% VAF
that were previously StrandBias-filtered now PASS when
`sv_strand_bias_threshold` is 1.0 (the default); SNV/indel strand bias
behaviour is unchanged; tests cover both paths.

## Steps

1. Read `kam-call/src/caller.rs` in full, focusing on `CallerConfig` and
   `assign_filter`.

2. Add the field to `CallerConfig`:
   ```rust
   /// Fisher p-value below which strand bias is flagged for SV types.
   ///
   /// Defaults to 1.0 (disabled). Inversion junction reads are structurally
   /// strand-biased due to directional path walking in the de Bruijn graph.
   /// The standard `strand_bias_threshold` is inappropriate for SV paths.
   ///
   /// Set to 0.0 to apply the same threshold as SNVs/indels.
   pub sv_strand_bias_threshold: f64,
   ```
   Default value: `1.0`.

3. In `assign_filter`, replace the strand bias check with:
   ```rust
   let eff_strand_bias_threshold = if is_sv {
       config.sv_strand_bias_threshold
   } else {
       config.strand_bias_threshold
   };
   if strand_bias_p < eff_strand_bias_threshold {
       return VariantFilter::StrandBias;
   }
   ```

4. Add or update tests:
   - An INV call with biased strands that passes when
     `sv_strand_bias_threshold = 1.0`.
   - The same call that is filtered when `sv_strand_bias_threshold = 0.01`.
   - Confirm an SNV call with biased strands is still filtered.

5. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- Setting the default to 1.0 means strand bias is never flagged for SV types.
  This is intentional: the strand imbalance is structural, not artefactual.
- The CLI flag to expose this value is in SV-STRAND-002.
