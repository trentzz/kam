# Review Cycle: Sensitivity Analysis
**Code version**: commit 87009ad
**Date**: 2026-03-20

---

## Code Review Findings

### Bug 1 (Critical): Strand bias test uses summed k-mer counts, not molecule counts

**Location**: `kam-pathfind/src/score.rs:161-162`, `kam-call/src/caller.rs:158-163`

`PathEvidence.total_simplex_fwd` and `total_simplex_rev` are computed as the **sum** of `n_simplex_fwd`/`n_simplex_rev` across all k-mers in the path. For a 100 bp target with k=31, the alt path has ~70 k-mers. Each alt molecule contributes to approximately 31 of the unique variant k-mers and ~39 shared k-mers. So with 5 alt forward-strand molecules:

- `total_simplex_fwd ≈ 5 × 70 = 350` (not 5)
- `total_simplex_rev = 0`

The Fisher exact test then sees (350, 0, ref_sum, ref_sum) instead of (5, 0, ref_mol, ref_mol). With inflated sample sizes, even a 5:0 forward:reverse split from a real variant with very few alt molecules is reported as extreme strand bias (p ≈ 0) and filtered.

A correct strand bias test requires per-molecule strand counts, not per-k-mer counts. The minimum of `n_simplex_fwd` across all path k-mers approximates the number of forward-strand molecules covering the whole path. Adding `min_simplex_fwd` and `min_simplex_rev` to `PathEvidence` and using these in the strand bias test will give properly calibrated p-values.

**Impact**: This bug causes many true positive calls to be filtered as `StrandBias`. Any real variant with few alt molecules (0.1–0.5% VAF) where all happen to be on one strand — common at low coverage — is incorrectly filtered. This is likely the dominant sensitivity loss factor.

### Bug 2 (Critical): Duplex orientation comparison fails for genuine duplex pairs

**Location**: `kam-assemble/src/assembler.rs:228-230`, `kam-assemble/src/fingerprint.rs:128-130`

For a genuine duplex pair, the forward template and reverse template produce fingerprints related by a 32-bit rotation:

```
fp_fwd = [enc(T[0:8]) | enc(T[-8:]) | enc(RC(T[-8:])) | enc(RC(T[0:8]))]
fp_rev = [enc(RC(T[-8:])) | enc(RC(T[0:8])) | enc(T[0:8]) | enc(T[-8:])]
       = rotate_left(fp_fwd, 32)
```

The current `fingerprints_compatible(fp1, fp2)` computes `(fp1 ^ fp2).count_ones() <= 4`. For a genuine duplex pair, `fp1 ^ fp2` has approximately 32 bits set (the upper and lower 32 bits encode different sequences), so the comparison always returns false. Forward and reverse templates of the same molecule are therefore split into separate fingerprint groups and never combined into a duplex molecule.

**Fix**: Also check `fingerprints_compatible(fp1, fp2.rotate_left(32))` in `split_by_fingerprint`. This allows genuine duplex pairs to be grouped together. Once grouped, `build_molecule` already has the logic to produce a duplex consensus from forward and reverse strand reads.

**Impact**: Once fixed, duplex rates should substantially increase at all depths. With more duplex evidence, more variants will have duplex support, improving confidence scores and the strand bias power.

### Bug 3 (Minor): Graph module doc comment says "canonical k-mers" but graph stores raw k-mers

**Location**: `kam-pathfind/src/graph.rs` doc comment, `kam-pathfind/src/score.rs:8-11`

The graph.rs module says "Nodes are canonical k-mers" but `score.rs` explicitly states nodes are raw k-mers and canonicalization happens at lookup. The score.rs description is correct. Fix the graph.rs doc comment.

### Bug 4 (Minor): Introduction section has PLACEHOLDER strings

**Location**: `docs/paper/sections/01_introduction.tex:28`

Line contains `\texttt{PLACEHOLDER_SENSITIVITY_LOW_VAF}` and similar. These must be filled with real numbers or removed before any further sharing of the paper.

---

## Paper Review Findings (Opus)

### Critical
1. PLACEHOLDER strings in Section 1 (introduction promises performance not yet achieved)
2. Strand bias methodology wrong in Section 4.6.3 (paper describes per-molecule counts; code uses per-k-mer sums)
3. No head-to-head comparison against any other tool

### Major
4. Abstract projects "0.1–0.5% VAF" range from a single-line duplex fix — speculative without data
5. "Precision above 0.55" framed as high — nearly half of calls are false positives
6. Graph node representation inconsistency between paper and graph.rs doc comment

### Minor
7. `k` used for both k-mer length and alt molecule count (Equation 8)
8. Empty appendix (`% TODO`)
9. ASCII-art pipeline figure in Section 2.3 not publication-quality
10. Duplex fix and O(n²) grouping listed as "future work" but are bugs/limitations
11. Reproducibility paragraph in Discussion reads as padding — move to Method
12. "single-line fix" in abstract is implementation detail, not abstract material

---

## Analysis and Plan

### Root cause of sensitivity ceiling at 56–59%

Two bugs compound to suppress sensitivity:

1. The strand bias overcounting makes the test ~70× more sensitive than intended. A 2:0 forward:reverse split among 2 alt molecules gives p ≈ 0.5 (not significant). The same split with summed k-mer counts gives p ≈ 10⁻⁵⁰ (filtered). Any low-VAF variant with few alt molecules is likely to appear strand-biased by chance, and this bug means it will be filtered even when the underlying molecule evidence is genuinely unbiased.

2. The duplex orientation bug prevents ~half the expected duplex molecules from being identified. Duplex evidence is both a direct confidence signal and the basis for the `LowDuplex` filter. With no duplex molecules, the `min_confidence = 0.99` threshold must be reached entirely from simplex evidence, which is harder at low alt counts.

### Changes planned

**Priority 1** — Fix strand bias overcounting:
- Add `min_simplex_fwd` and `min_simplex_rev` to `PathEvidence` in `score.rs`
- Use these in `call_variant` for the strand bias test
- Update all tests and doc examples

**Priority 2** — Fix duplex orientation:
- Add `fingerprints_duplex_compatible(fp1, fp2)` to `fingerprint.rs` that also checks the 32-bit rotation
- Use it in `split_by_fingerprint` in `assembler.rs`
- Add tests for genuine duplex pair grouping

**Priority 3** — Paper corrections:
- Fill or remove introduction placeholders
- Correct the strand bias methodology description
- Remove "single-line fix" from abstract
- Move duplex bug and O(n²) grouping from future work to evaluation/limitations
- Fix variable collision (`k` for k-mer length vs alt count)
