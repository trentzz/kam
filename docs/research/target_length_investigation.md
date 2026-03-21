# Investigation: Target Window Length and Sensitivity

**Date**: 2026-03-22
**Status**: Benchmark running. 50bp and 60bp results corrected after anchor overlap guard bug fix.

---

## Motivation

The primary cause of sensitivity loss at 2% VAF is end-anchor coverage failure for indels and
per-target coverage imbalance for SNVs (see `sensitivity_investigation.md`). Both causes worsen
when the target window is short (fewer reads span the anchor) and improve when the window is
longer (more reads span the anchor, reducing the fraction of coverage-gap failures).

The panel uses 100bp target windows. This investigation measures whether extending to 120bp,
150bp, or 200bp would meaningfully improve sensitivity, and whether shorter windows would
degrade it.

---

## Anchor Overlap Guard Bug (50bp and 60bp showed sens=0.000)

### Symptom

Initial benchmark showed `sens=0.000` for all 50bp and 60bp samples despite
`start_missing=0, end_missing=0` (both anchors present in the graph).

### Root Cause

`run.rs` contained an anchor overlap guard:

```rust
let end_anchor_pos = target_seq.len().saturating_sub(k + end_offset);
if start_offset + k > end_anchor_pos {
    n_walk_no_paths += 1;
    ...
    continue;
}
```

For a target of length `L` with kmer size `k=31` and both offsets = 0:
- `end_anchor_pos = L - 31`
- Guard fires when: `0 + 31 > L - 31`, i.e. when `L < 62`

This unconditionally skips all targets shorter than `2k = 62bp`. For 50bp targets,
all 375 targets are skipped → `no_paths=375` → `sens=0.000`.

### Fix

The guard intended to catch a degenerate case where aggressive soft anchoring from both
ends caused the effective start and end anchors to be the same k-mer. The correct check is
whether the two anchor k-mers are identical (`start_raw == end_raw`), not whether their
genomic windows overlap. Overlapping windows are fine: the two k-mers are still distinct
and the DFS walk finds valid paths between them.

The guard was replaced with:

```rust
if start_raw == end_raw {
    // Anchors are the same k-mer; no variant path can be distinguished.
    n_walk_no_paths += 1;
    ...
    continue;
}
```

This allows short targets (≥ k+1 = 32bp) while still skipping the pathological case.

### Lesson

When a diagnostic shows `start_missing=0, end_missing=0, no_paths=375`, the walk is
being skipped BEFORE it starts, not failing during traversal. Always check for early-exit
guards in the path between anchor validation and `walk_paths()`.

---

## Results (2M reads, --max-vaf 0.35 --min-family-size 2)

Results focus on the 15ng 2% VAF and 30ng 2% VAF samples as the primary sensitivity
benchmark. SNV and indel sensitivity are reported separately because they have different
sensitivity ceilings and respond differently to target length changes.

*(Table to be filled once all benchmark runs complete.)*

| Target length | 15ng 2% sens | 15ng 2% SNV | 15ng 2% indel | 30ng 2% sens | Notes |
|--------------|-------------|-------------|---------------|-------------|-------|
| 50bp  | TBD | TBD | TBD | TBD | Anchor overlap guard bug fixed |
| 60bp  | TBD | TBD | TBD | TBD | Anchor overlap guard bug fixed |
| 70bp  | TBD | TBD | TBD | TBD | |
| 80bp  | TBD | TBD | TBD | TBD | |
| 90bp  | TBD | TBD | TBD | TBD | |
| 100bp | 0.613 | 0.800 | 0.388 | 0.592 | Baseline (v10) |
| 120bp | TBD | TBD | TBD | TBD | |
| 150bp | TBD | TBD | TBD | TBD | |
| 200bp | TBD | TBD | TBD | TBD | |

---

## Minimum Target Length for Anchored DFS

For the anchored DFS walk to detect a variant, both the start and end anchor k-mers must be
shared between the reference and variant paths. This requires the variant position to fall
outside the first k and last k bases of the target sequence.

For a target of length L (note: FASTA files are L+1 bp due to inclusive coordinate generation)
with a variant centred at position L/2 and k=31:

- The start anchor k-mer covers positions 0..31. The variant must be at position ≥ 31 to
  avoid affecting the start anchor.
- The end anchor k-mer covers positions L-31..L. The variant must be at position ≤ L-32 to
  avoid affecting the end anchor.
- The detectable range is positions 31 to L-32.
- This range is non-empty when L > 62 (actual sequence length: 63bp).

For the actual benchmarked files (all off by +1 due to coordinate convention):

| Target label | Actual sequence length | Min detectable position | Max detectable position | Centred variant at |
|-------------|----------------------|------------------------|------------------------|---------------------|
| 50bp  | 51 bp | 31 | 20 (empty!) | 25 (not detectable) |
| 60bp  | 61 bp | 31 | 30 (empty!) | 30 (not detectable) |
| 70bp  | 71 bp | 31 | 40 | 35 (detectable, margin=5) |
| 80bp  | 81 bp | 31 | 50 | 40 (detectable, margin=10) |
| 90bp  | 91 bp | 31 | 60 | 45 (detectable, margin=15) |
| 100bp | 101 bp | 31 | 70 | 50 (detectable, margin=20) |
| 120bp | 121 bp | 31 | 90 | 60 (detectable, margin=30) |
| 150bp | 151 bp | 31 | 120 | 75 (detectable, margin=45) |
| 200bp | 201 bp | 31 | 170 | 100 (detectable, margin=70) |

The 50bp and 60bp targets are fundamentally undetectable with k=31. Even after the anchor
overlap guard bug was fixed, these benchmarks correctly return sens=0.000. This is not a
bug — it is a genuine algorithmic limitation of anchored DFS when target_len < 2k.

The minimum usable target length for k=31 is 63bp (the variant must be at least 32 bases
from each anchor). In practice, some targets will have off-centre variants, and coverage
imbalance further reduces sensitivity, so effective minimum is larger.

---

## Expected Effects of Target Length

### Why longer targets should improve sensitivity

1. **End-anchor coverage**: the end anchor k-mer sits at the last k=31 bases of the target.
   Longer targets push the end anchor further from the panel probe boundaries. Reads that
   currently cut off before the 100bp end anchor would span the anchor in a 150bp window,
   increasing the fraction of targets with a valid end anchor.

2. **Indel walk completion**: for indels of length L, the reference path is L k-mers shorter
   (deletion) or longer (insertion) than expected. The walk needs reads spanning from the
   variant junction to the end anchor. A longer target shifts the end anchor further right,
   reducing the fraction of indel paths truncated before reaching the anchor.

3. **Coverage imbalance at boundaries**: targets with low coverage at their exact 100bp
   boundaries would gain extra covered k-mers if the window were wider.

### Why very short targets hurt sensitivity

For short targets (< 62bp for k=31), the anchor windows overlap, but the walk itself is
still valid (as shown by the bug investigation above). The real problem with short targets
is that a variant at the centre of a 50bp window is only ~20bp from either anchor, and
reads providing coverage of that variant must also span to the anchor — a shorter span
than in a 100bp window. This reduces the fraction of reads that provide useful variant
evidence at the walk boundaries.

---

## Key Measurement: Indel Sensitivity vs Target Length

Indel sensitivity is expected to be the metric most responsive to target length. The 61–68%
indel false-negative rate at 100bp is dominated by end-anchor displacement. Longer windows
reduce displacement.

A clean result would show indel sensitivity rising monotonically with target length up to
some plateau, while SNV sensitivity shows a smaller or flat response (SNVs do not have
end-anchor displacement; their missed calls come from per-target coverage imbalance which
is less target-length-dependent).

---

## Reference

- `sensitivity_investigation.md` — primary sensitivity ceiling analysis
- `anchor_missing_investigation.md` — soft anchor fix and its limited impact
- `hash_partition_umi_grouping.md` — Phase 1 and 2 UMI grouping improvements
