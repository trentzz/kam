# SV-EXP-006: Fusion Classification Logic

**Epic**: SV-EXPAND (docs/claudetracking/overallplans/SV-EXPAND.md)
**Priority**: high
**Depends on**: sv_exp_005_fusion_junctions.md
**Status**: todo

## Goal

Implement the classifier logic that identifies a walked path as a `Fusion` when
its k-mers span two distinct `target_id` values. Done looks like: a path built
from junction k-mers spanning partner A and partner B is classified as
`Fusion`; a normal intra-target path is not affected.

## Steps

1. Read `kam-call/src/classify.rs` and the fusion detection design document
   (`docs/research/fusion_detection_design.md`).

2. In the classification step, after the path is walked:
   - Collect the set of `target_id` values from all k-mers in the path.
   - If the set contains exactly two distinct `target_id` values and the path
     includes at least one junction k-mer (annotated in SV-EXP-005), classify
     as `Fusion`.

3. Add a `Fusion` arm to all relevant match blocks in `classify.rs`.

4. Record the two partner `target_id` values in the variant record for use in
   VCF output (SV-EXP-009).

5. Write unit tests:
   - Path with k-mers from two targets → `Fusion`.
   - Path with k-mers from one target → not `Fusion`.
   - Path with k-mers from three targets (ambiguous) → `Fusion` or discarded;
     choose and document the behaviour.

6. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- The implementation must follow the design from SV-EXP-004 exactly. If the
  design changes during implementation, update the design document first.
- Do not implement BND VCF output here; that is SV-EXP-009.
