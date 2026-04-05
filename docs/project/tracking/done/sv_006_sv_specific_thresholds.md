# SV-006: SV-Type-Specific Caller Thresholds

## Problem

The current caller applies uniform `min_confidence = 0.99` and
`min_alt_molecules = 2` to all variant types. This is too strict for
structural variants.

A 100 bp deletion or inversion with 2 supporting molecules has near-zero
probability of arising from a single-position background error. The binomial
confidence of 0.981 with 2 alt molecules at the DUP junction correctly
reflects strong signal; rejecting it at the 0.99 threshold is overly
conservative.

In the SV benchmark (sv_sensitivity_report.md), TandemDuplication fails at
0.5% VAF with 2 molecules and confidence 0.981. Only the threshold is
preventing a PASS call.

## Fix

Add SV-specific thresholds to `CallerConfig`:

- `sv_min_confidence: f64` — confidence threshold for SV types (LargeDeletion,
  TandemDuplication, Inversion). Default: 0.95.
- `sv_min_alt_molecules: u32` — minimum alt molecule count for SV types.
  Default: 1.

Expose both via CLI flags: `--sv-min-confidence` and `--sv-min-alt-molecules`.

In `assign_filter`, accept the `VariantType` and switch to SV-specific
thresholds when the type is an SV.

## Test Plan

1. Extend the varforge SV benchmark with a 0.25% VAF level.
2. Re-run all VAF levels (0.25%, 0.5%, 1%, 2%, 5%).
3. Show DUP at 0.5% VAF now PASS (confidence 0.981 > 0.95).
4. Confirm no new false positives in discovery mode.
5. Confirm monitoring mode still eliminates all FPs.

## Expected Result

- DUP PASS at 0.5% VAF (previously LowConfidence).
- DUP LowConfidence at 0.25% VAF (1 molecule, confidence ~0.83 < 0.95).
- All other SV types unchanged (already PASS at all tested VAFs).
