# Sensitivity Investigation: Why kam Misses Variants

**Date**: 2026-03-21
**Benchmark**: v7 baseline (2M reads, `--max-vaf 0.35 --min-family-size 2 --target-variants`)
**Sensitivity at 2% VAF**: 51.7%–61.3%

---

## Panel and depth characteristics

- 375 truth variants (205 SNV, 170 indel) across 375 target regions (100bp each)
- At 2M reads: 293K–416K assembled molecules per sample
- Molecules per target: ~780–1110 (varies by concentration)
- Mean alt molecules per target at 2% VAF: ~16–22 (Poisson mean = 0.02 × molecules_per_target)
- Mean alt molecules per target at 0.5% VAF: ~4–6

---

## Hypothesis 1: min_alt_molecules filter (disproved)

Initial hypothesis: setting `min_alt_molecules=2` filters out legitimate 1-molecule calls,
which account for ~22% of targets at 2% VAF under a Poisson model.

Test: benchmark v8 (`--min-alt-molecules 1`, default confidence=0.99) and v9
(`--min-alt-molecules 1 --min-confidence 0.95`).

Results at 2% VAF:
- 15ng: v7=0.613, v9=0.613 (no change)
- 30ng: v7=0.592, v9=0.597 (+0.5pp)
- 5ng: not yet measured (expected similar)

Conclusion: parameter tuning yields ≤0.5pp improvement. The filters are not the
primary bottleneck.

**Why the hypothesis was wrong**: at 2% VAF with ~780 molecules per target,
the Poisson mean is ~16 alt molecules per target. P(k=0) ≈ 10^-7. P(k=1) ≈ 10^-6.
Almost every truth variant has ≥10 alt molecules. The filter at k≥2 rarely fires.

---

## Actual bottleneck: graph walk failures

The sensitivity ceiling is set by how many truth variant paths the DFS walk finds.
At 2% VAF, SNV sensitivity is 69–80% and indel sensitivity is 31–39%.

### SNV walk failures (20–31% of SNVs missed)

Most missed SNVs arise from two causes:

1. **Per-target coverage imbalance**: molecule depth is not uniform across targets.
   GC content and PCR amplification bias cause some targets to have far fewer
   molecules than the panel average. A target with 5 molecules at 2% VAF has only
   ~0.1 alt molecules on average; most runs will find zero, and the walk either
   finds no alt path or fails the confidence filter.

2. **Anchor coverage gaps**: the DFS walk requires the start AND end anchor k-mers
   to be present in the graph. For some targets, the reference end-anchor k-mer
   may be absent (reads cut off before the target end due to short inserts).

### Indel walk failures (61–68% of indels missed)

The dominant cause: end-anchor displacement.

For an indel of length L, the variant path is L k-mers longer (insertion) or
shorter (deletion) than the reference path. The reference end-anchor k-mer
(last k-mer of the 100bp target) sits at position 70 - L relative to the
post-indel sequence (for deletion) or 70 + L (for insertion).

For the walk to complete, reads must span from the variant junction to the
reference end-anchor. This is the "anchor coverage" requirement. At lower VAF,
variant reads may not cover the full span from junction to end-anchor. For
larger indels (L > 10), the required span is long enough that many reads cut
off before reaching the anchor.

This is a structural limitation of the current anchored-DFS walk design.

---

## What would actually improve sensitivity

### Short-term (incremental, no algorithm change)

None. Parameter tuning yields at most 1–2pp improvement.

### Medium-term (algorithm change, non-trivial)

1. **Relaxed end-anchor matching for indels**: instead of requiring the exact last
   k-mer of the target, allow the walk to terminate within a window of ±L k-mers
   of the end anchor. This would recover indel paths where reads stop just short
   of the end anchor.

2. **Read-level k-mer inclusion**: the current graph uses MOLECULE consensus reads.
   Including raw reads (before consensus) would provide more k-mer coverage at
   the cost of higher noise. With min_family_size=2, most consensus reads already
   incorporate ≥2 reads, so the coverage gain would be marginal.

3. **Target extension**: extend target windows from 100bp to 150bp. This shifts
   the end anchor further from most variants, reducing the fraction of cases where
   reads cut off before the end anchor.

### Long-term (depth increase)

The O(n²) UMI grouping bottleneck limits practical depth to ~2M read pairs.
Hash-partition grouping (O(n)) would allow 10–20M read pairs, giving 5–10x more
alt molecules per target and recovering many of the coverage-limited FNs.

---

## Summary

The 52–61% sensitivity at 2% VAF is set by walk failures, primarily for indels
(61–68% missed) and some SNVs (20–31% missed). These failures are caused by
end-anchor displacement (indels) and per-target coverage imbalance (SNVs), not
by the confidence or alt molecule filters.

Meaningful sensitivity improvement requires either:
- Algorithm changes to the end-anchor walk (medium-term)
- Deeper sequencing via hash-partition UMI grouping (long-term)

The current v7 parameter set (`--max-vaf 0.35 --min-family-size 2
--target-variants`) is already near the sensitivity ceiling for this algorithm
and read depth.
