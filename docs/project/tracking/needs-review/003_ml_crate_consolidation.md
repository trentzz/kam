# Consolidate kam-ml into kam-call?

**Category**: architecture
**Related epic**: QUALITY

## Context

`kam-ml` is 224 lines with one public type (`MlScorer`) and one private helper (`variant_type_str` — which duplicates the existing `Display for VariantType` in kam-call). The crate was created per `docs/project/devmanual/ml_integration_design.md` to separate ONNX runtime from core types. In practice:

- kam-ml depends on kam-call (for `VariantCall`/`VariantType`).
- kam-call depends on kam-pathfind for `PathEvidence`.
- `kam-ml` is gated behind `feature = "ml"` already, so binary size concerns are already addressed at the crate level.

Arguments for consolidation:
- Removes a workspace member.
- Kills `variant_type_str` duplication.
- Single import surface for downstream users.

Arguments against:
- Clean boundary lets you swap ONNX for a different runtime.
- Future ML features (training-time utilities, calibration) would bloat kam-call.

## Options

1. **Merge into `kam-call::ml` module, feature-gated**. Deletes kam-ml.
2. **Keep separate**. Dedup `variant_type_str` via a shared helper in kam-call.
3. **Grow kam-ml**. Move feature extraction (currently inlined in the scorer) and calibration utilities into kam-ml, so the crate earns its weight.

## Recommendation

Option 2 for now; revisit after v2 of the ML model ships. The crate boundary is cheap and paying off (clean feature gating). Fix the duplication by adding a `kam_call::caller::variant_type_tag() -> &'static str` helper and re-exporting.
