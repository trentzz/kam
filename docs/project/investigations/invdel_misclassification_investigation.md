# InvDel Misclassification as Deletion

**Date**: 2026-03-27
**Status**: Complete.

---

## Symptom

InvDel variants were correctly detected by the pathfinder (alt_found=1, correct junction path recovered) but appeared in benchmark output as `Deletion` with filter `StrandBias`. No InvDel calls reached PASS or LowConfidence at any VAF level.

## Hypothesis

`classify_variant` in `kam-call/src/caller.rs` gates InvDel detection behind a net-length-change check. An InvDel with a large deleted segment and a shorter inverted-inserted segment can have a small net length change. If the net change falls below `SV_LENGTH_THRESHOLD` (50 bp), the function returns `Deletion` without ever testing for the inversion pattern. This would then trigger `assign_filter` with SNV/indel strand bias thresholds (0.01) rather than the SV threshold (1.0, effectively disabled), causing StrandBias to filter real InvDel calls.

## Root Cause

Confirmed. The synthetic benchmark uses an InvDel with 79 bp deleted and 60 bp inverted-inserted. The net length change is 79 - 60 = 19 bp. The check `if del_len >= SV_LENGTH_THRESHOLD` evaluated to `19 >= 50 = false`, so the function returned `Deletion` without calling `alt_seq_has_inversion_relative_to_ref`.

`assign_filter` then received `VariantType::Deletion` and applied the SNV/indel strand bias threshold (0.01). The real InvDel signal had imperfect strand balance (as expected for a low-VAF SV), so it was filtered as StrandBias.

The inversion detection function itself is correct. It internally requires both the ref central segment and the alt central segment to be at least `SV_LENGTH_THRESHOLD` (50 bp) before accepting an inversion call, so it would not falsely classify a short deletion as InvDel. The problem was that the function was never called when the net length change was small.

## Fix

Moved the `alt_seq_has_inversion_relative_to_ref` call before the `del_len >= SV_LENGTH_THRESHOLD` gate in `classify_variant`. If the inversion test passes, the function returns `InvDel` immediately, regardless of net length change. Net-length gating for the `LargeDeletion` vs `Deletion` distinction is preserved for all non-InvDel paths. Commit `0a9406f`.

## Measured Result

After the fix, InvDel variants are correctly classified and reach PASS with SV-specific thresholds:

| VAF  | Before fix         | After fix    |
|------|--------------------|--------------|
| 0.5% | Deletion/StrandBias | **PASS**     |
| 1%   | Deletion/StrandBias | **PASS**     |
| 2%   | Deletion/StrandBias | **PASS**     |
| 5%   | Deletion/StrandBias | **PASS**     |

No regressions in Deletion or LargeDeletion classification. The gate order change only affects paths where `alt_seq_has_inversion_relative_to_ref` returns true, which requires both central segments to be at least 50 bp.
