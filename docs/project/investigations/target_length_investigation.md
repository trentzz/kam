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

Note: target FASTA files contain sequences that are L+1 bp long (e.g., "50bp" files contain
51bp sequences) due to inclusive coordinate generation. The actual sequence length is used in
the minimum-length analysis above.

| Target label | Actual seq len | 15ng 2% sens | 15ng 2% SNV | 15ng 2% indel | 30ng 2% sens | 5ng 2% sens | Pathfind ms | Notes |
|-------------|---------------|-------------|-------------|---------------|-------------|-------------|-------------|-------|
| 50bp  | 51 bp  | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | ~590  | Below 2k minimum; algorithmic limit |
| 60bp  | 61 bp  | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | ~600  | Below 2k minimum; algorithmic limit |
| 70bp  | 71 bp  | 0.603 | 0.805 | 0.359 | 0.581 | 0.507 | 940   | Just above minimum; narrow margin |
| 80bp  | 81 bp  | 0.616 | 0.805 | 0.388 | 0.595 | 0.517 | 1224  | **Peak sensitivity** |
| 90bp  | 91 bp  | 0.613 | 0.800 | 0.388 | 0.592 | 0.517 | 1326  | Same as 100bp |
| 100bp | 101 bp | 0.613 | 0.800 | 0.388 | 0.592 | 0.517 | 1726  | Baseline (current panel) |
| 120bp | 121 bp | 0.600 | 0.780 | 0.382 | 0.581 | 0.512 | 3132  | Spurious paths; slightly worse |
| 150bp | 151 bp | 0.584 | 0.761 | 0.371 | 0.568 | 0.507 | 9165  | Notably worse; 5.3× pathfind time |
| 200bp | 201 bp | 0.483 | 0.620 | 0.318 | 0.531 | 0.499 | 27103 | Severely degraded; 15.7× pathfind time |

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

## Root Cause: Why Longer Targets Degrade Sensitivity

For 120bp and beyond, sensitivity declines despite reads easily spanning the full target.
The cause is graph complexity and the `max_paths=100` walk limit.

**Pathfind diagnostics (15ng 2% VAF):**

| Target | no_paths | ref_only | alt_found | time_ms |
|--------|---------|----------|-----------|---------|
| 100bp  | 67      | 35       | 273       | 1701    |
| 120bp  | 69      | 29       | 277       | 3228    |
| 150bp  | 75      | 22       | 278       | 9252    |
| 200bp  | 88      | 11       | 276       | 26776   |

`alt_found` stays high (and even increases) with longer targets — the walk succeeds in
finding alt paths for more targets. But sensitivity falls because the paths found do not
match the truth. Longer targets have more k-mers in the graph, more branching points, and
more spurious variant paths. Once the DFS fills the `max_paths=100` budget with spurious
paths, the true variant path is not reached, and the call scores as non-truth.

The truth-matching call rate falls from 230/273 (84%) at 100bp to 225/277 (81%) at 120bp
to approximately 219/278 (79%) at 150bp. Each alt path found has a lower probability of
being the truth variant.

The `no_paths` increase at 200bp (88 vs 67 at 100bp) reflects a different mechanism:
anchor coverage failures. For 201bp targets with 150bp reads, the anchors at positions
0–31 and 170–200 cannot be fully covered by a single read. Coverage at the far boundaries
degrades, causing the walk to fail before starting.

**Practical impact:** the 200bp pathfind time of 27 seconds per sample is 15.7× the
1.7-second 100bp pathfind. Even without the sensitivity penalty, 200bp targets are
unusable for high-throughput pipelines.

---

## Key Measurements

**Indel sensitivity vs target length**: Indel sensitivity rises from 0.359 at 70bp to 0.388
at 80bp (same as 100bp), then declines for longer targets (0.371 at 150bp, 0.318 at 200bp).
The improvement from 70bp to 80bp confirms the theoretical prediction that a wider window
gives more coverage for shifted end-anchor positions. The plateau at 80–100bp and decline
beyond indicates that the spurious path and coverage gap effects dominate.

**SNV sensitivity vs target length**: SNV sensitivity is roughly constant at 0.800–0.805
from 70bp to 100bp, then falls to 0.761 at 150bp and 0.620 at 200bp. SNVs do not have
end-anchor displacement (they don't change path length), so the gain from 70bp to 80bp is
smaller than for indels. The decline at 120bp+ is entirely due to spurious path crowding.

**Conclusion**: the optimal target length for k=31 with 150bp paired-end reads is
**80–100bp**. The current panel uses 100bp, which is near-optimal.

Recommended setting: do NOT change to longer targets. The current 100bp panel design is
within the plateau and avoids the spurious path problem. If re-designing for a future panel,
80bp would give marginally better indel sensitivity with marginally faster pathfind.

---

## Reference

- `sensitivity_investigation.md` — primary sensitivity ceiling analysis
- `anchor_missing_investigation.md` — soft anchor fix and its limited impact
- `hash_partition_umi_grouping.md` — Phase 1 and 2 UMI grouping improvements
