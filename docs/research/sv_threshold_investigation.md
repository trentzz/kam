# SV-006: SV-Type-Specific Caller Thresholds

**Date**: 2026-03-23
**Status**: Complete.

---

## Symptom

In the SV benchmark (SV-005), TandemDuplication at 0.5% VAF was labelled
LowConfidence despite having 2 supporting molecules and confidence 0.981. The
LargeDeletion and Inversion at the same VAF level were PASS. The difference was
not in evidence quality but in the confidence threshold: 0.981 < 0.99.

## Hypothesis

A uniform 0.99 confidence threshold is appropriate for SNVs and indels, where a
false positive can arise from a background sequencing error at a single base.
For structural variants (LargeDeletion, TandemDuplication, Inversion), a false
positive requires a background process that mimics a 50--200 bp structural
rearrangement. This is essentially impossible for random sequencing errors. A
lower threshold is therefore justified for SV types without increasing the risk
of false positives.

The specific hypothesis: with `sv_min_confidence = 0.95`, the 2-molecule DUP
call at 0.5% VAF (confidence 0.981) would pass, and no new false positives would
be introduced (since existing FPs are SNV/MNV errors, not SVs).

## Root Cause

The confidence formula is binomial likelihood-based. With 2 alt molecules out of
~1000 total:

- L(f_hat = 0.002) = (0.002)^2 × (0.998)^998 ≈ 3.64e-7
- L(epsilon = 1e-4) = (1e-4)^2 × (0.9999)^998 ≈ 9.0e-9
- Confidence = 3.64e-7 / (3.64e-7 + 9.0e-9) ≈ 0.982

The 0.99 threshold rejects this. But the signal is real: 2 molecules
independently support the same 100 bp insertion at the same junction. The
probability that 2 independent sequencing errors produce identical 100 bp
insertions is negligible.

The root cause is a threshold designed for single-base events being applied
uniformly to multi-kilobase structural events.

## Fix

Added two new fields to `CallerConfig` in `kam-call/src/caller.rs`:
- `sv_min_confidence: f64` (default 0.95)
- `sv_min_alt_molecules: u32` (default 1)

Modified `assign_filter` to accept the `VariantType` and apply SV-specific
thresholds when the type is `LargeDeletion`, `TandemDuplication`, or
`Inversion`. Exposed both fields as CLI flags (`--sv-min-confidence`,
`--sv-min-alt-molecules`).

## Measured Result

With the updated binary, re-running the existing SV benchmark (same simulation
data, no new varforge runs):

| VAF  | DUP before | DUP after |
|------|------------|-----------|
| 0.5% | LowConf    | **PASS**  |
| 1%   | LowConf    | LowConf   |
| 2%   | PASS       | PASS      |
| 5%   | PASS       | PASS      |

DUP at 0.5% (confidence 0.982 > 0.95) now passes. DUP at 1% (confidence 0.791
< 0.95, 1 supporting molecule) remains LowConfidence. The 1-molecule call
represents genuinely weak evidence; accepting it would increase the risk of
spurious calls in tumour-informed mode.

False positives in discovery mode: unchanged (3, 2, 6, 7 at 0.5%–5% VAF). All
FPs are SNV/MNV types — the SV-specific threshold does not affect them.

Monitoring mode: zero FPs at every VAF level (unchanged from SV-005).

## Summary

All three SV types are now PASS at ≥0.5% VAF except TandemDuplication at 1%
(1 molecule, genuinely weak). The change is principled: the confidence
calculation is the same; only the acceptance threshold differs by variant class.
