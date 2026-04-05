# Fusion Detection VAF=0 Root Cause

**Date**: 2026-03-27
**Status**: Complete.

---

## Symptom

Fusion detection returned VAF=0.0 and molecules=0 at every VAF level in the benchmark. The pathfinder was finding fusion junction paths (alt_found=1), but the caller reported zero supporting molecules regardless of input VAF (0.5%–5%).

## Hypothesis

Two independent hypotheses were tested.

**Hypothesis 1**: `score_and_rank_paths` in `kam-pathfind/src/score.rs` was skipping the reference path during molecule evidence computation. For fusion calling, the fusion junction path is matched against the fusion target sequence and marked `is_reference=true`. If the evidence loop skips that path, `mean_variant_specific_molecules` stays at 0.0, and the caller sees k=0.

**Hypothesis 2**: `find_partner_depth` in `kam/src/commands/run.rs` was matching fusion partners by chromosome name substring. When both fusion partners lie on the same chromosome, HashMap iteration order is non-deterministic. Both A and B depth lookups could return the same entry, making the VAF denominator wrong by up to 2x. This alone would not cause VAF=0, but it would skew the result.

## Root Cause

Both hypotheses were confirmed.

**Bug 1 (critical)** in `kam-pathfind/src/score.rs`: `score_and_rank_paths` contained an early `continue` when `sp.is_reference` was true. This guard was intended to skip populating `mean_variant_specific_molecules` for normal reference paths, since those k-mers are not variant-specific. For fusion calling, the fusion junction path is the wanted signal. It matches the fusion target sequence, so it is marked `is_reference=true`. The early `continue` caused it to leave `vs_evidences` empty. The existing fallback sets `mean_variant_specific_molecules = mean_molecules` when `vs_evidences` is empty, but that fallback was never reached. Instead, `mean_variant_specific_molecules` was left at its default of 0.0. `call_fusion` in `fusion.rs` reads `ref_evidence.mean_variant_specific_molecules` for the molecule count, so k=0 and VAF=0 in every case.

**Bug 2 (secondary)** in `kam/src/commands/run.rs`: `find_partner_depth` iterated a HashMap and returned the first entry whose chromosome name contained the query string as a substring. For same-chromosome fusions, both the A and B locus queries matched the same chromosome string, so iteration could return the same entry for both. This made the depth denominator unreliable.

## Fix

**Bug 1**: Removed the `if sp.is_reference { continue; }` guard in `score_and_rank_paths`. The fallback logic (`mean_variant_specific_molecules = mean_molecules` when `vs_evidences` is empty) is correct and handles normal reference paths that have no variant-specific k-mers. Commit `70e18b6`.

**Bug 2**: Changed `find_partner_depth` to accept a locus midpoint coordinate alongside the chromosome name. The function now filters by exact chromosome match, then selects the target whose encoded coordinate range midpoint is nearest to the query midpoint. This correctly disambiguates two fusion partners on the same chromosome. Commit `e6c0e3e`.

## Measured Result

After both fixes, fusion calls match expected VAF at each level:

| VAF  | Before fix | After fix |
|------|------------|-----------|
| 0.5% | 0.0 (0 mol) | ~0.005 |
| 1%   | 0.0 (0 mol) | ~0.010 |
| 2%   | 0.0 (0 mol) | ~0.020 |
| 5%   | 0.0 (0 mol) | ~0.050 |

Bug 1 was the primary cause of VAF=0. Bug 2 affected the denominator but would not have produced zero on its own. Fixing Bug 1 alone restored correct molecule counts. Fixing Bug 2 removed the non-deterministic partner depth selection.
