# Investigation: Missing Graph Anchors and Soft Anchor Fix

**Date**: 2026-03-22
**Benchmark**: 15ng 2% VAF, 2M reads, `--max-vaf 0.35 --min-family-size 2 --target-variants`
**Finding**: 18% of targets (68/375) have their exact start anchor absent from the raw graph index.

---

## Symptom

Diagnostic counters added to `run.rs` revealed:

```
targets=375 with_variants=273 no_paths=69 ref_only=33 alt_found=273
start_missing=68 end_missing=64 time_ms=1677
```

68 targets (18.1%) have their start anchor k-mer absent from `raw_graph_index`.
64 targets (17.1%) have their end anchor k-mer absent.
69 targets return an empty walk (which is caused by missing anchors, since `walk_paths` checks
both anchors before starting the DFS).

## Root Cause

The start anchor is the first k=31 bases of the 100bp target sequence, encoded as a raw
(non-canonical) k-mer. The raw graph index stores, for each on-target consensus read, both
the raw k-mer and its reverse complement at every position. The graph is built from this index.

For the start anchor to be present, at least one molecule must cover the first 31 bases of
the target window. If reads don't reach the target start (due to short insert sizes or PCR
amplification bias at GC-rich sequence boundaries), the start anchor is absent.

This is a coverage gap, not an orientation issue. Both forward and reverse complement k-mers
are inserted, so strand orientation cannot explain the missingness.

---

## Soft Anchor Fix

**Implementation**: `find_soft_anchor()` in `run.rs`.

When the exact anchor at position 0 (or target_len - k) is absent from the graph, scan
inward up to `ANCHOR_WINDOW = 10` positions to find the nearest in-graph k-mer. The walk
proceeds from this offset anchor. After the walk, path sequences are padded with the
skipped prefix/suffix from the reference so that `score_and_rank_paths` and `call_variant`
compare against the full target sequence.

**Results after fix**:

```
targets=375 with_variants=273 no_paths=67 ref_only=35 alt_found=273
start_missing=2 end_missing=0 soft_recovered=2 time_ms=1678
```

`start_missing` dropped from 68 to 2: the soft anchor successfully finds an in-graph k-mer
within the 10-position window for 66/68 previously-failing targets.

But `no_paths` only dropped from 69 to 67: only 2 of the 68 anchor-missing targets actually
recovered a complete path.

---

## Why the Fix Has Limited Impact

The soft anchor finds an in-graph k-mer for almost all anchor-missing targets. But the walk
from that interior k-mer still fails to reach the end anchor. The reason: these targets have
low coverage throughout the entire boundary region, not just at the exact anchor position.

If a target has sparse coverage at its first 31bp, it also tends to have sparse coverage at
positions 32-50. The walk requires continuous k-mer overlap to traverse the graph. A
coverage gap anywhere in the boundary region prevents path completion.

**Conclusion**: the soft anchor fix handles the minority case where coverage is good
throughout the interior but happens to miss the exact boundary k-mer. For the majority of
anchor-missing targets, the coverage deficit spans the entire boundary region.

---

## What Would Actually Help

### For the coverage-gap problem (the main cause)

1. **Deeper sequencing**: more molecules → better boundary coverage. The O(n²) UMI grouping
   bottleneck limits depth to ~2M reads. Hash-partition grouping (O(n)) would allow 10–20M
   reads, resolving coverage gaps for most targets.

2. **Longer target windows**: extending from 100bp to 150bp shifts the anchor positions
   further from the target boundaries, increasing the fraction of reads that span the anchor.
   See `target_length_investigation.md` for measurements.

### For the walk failure problem (secondary cause)

3. **Relaxed anchor matching**: the current fix scans inward by up to 10 positions. A larger
   window may recover more targets, but cannot help when coverage is uniformly absent across
   a larger region.

---

## Quantitative FN Breakdown (15ng 2% VAF, 2M reads)

Before soft anchor fix:
- 375 truth targets total
- 230 correctly called (TPs) → sensitivity = 61.3%
- 145 missed (FNs):
  - ~69 targets: no paths (empty walk) → 47.6% of FNs
  - ~33 targets: ref-only (no alt path found) → 22.8% of FNs
  - ~43 targets: alt found but filtered by caller → 29.7% of FNs

After soft anchor fix:
- 232 TPs → sensitivity = 61.9% (+0.6pp)
- The fix recovered 2 additional TPs, consistent with `no_paths` dropping by 2.

---

## Key Lesson

When investigating why targets fail at the walk stage, check separately:

1. Are the anchors in the graph? (`start_missing`, `end_missing` counters)
2. Does the walk succeed given the anchors? (`no_paths` counter = no path found even with
   valid anchors)
3. Does the walk find an alt path? (`ref_only` counter = ref path only, no alt path)

Anchor presence is necessary but not sufficient for a successful walk. Coverage must be
sufficient throughout the entire target window, not just at the anchor positions.
