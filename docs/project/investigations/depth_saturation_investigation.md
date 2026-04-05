# Investigation: Read Depth vs Sensitivity — Saturation at 2M Reads

**Date**: 2026-03-22
**Status**: Complete

---

## Symptom

A read depth sweep from 200K to 10M reads (15ng, 2% VAF, k=31, 100bp targets) showed
sensitivity plateauing sharply at 2M reads and not improving at 5M or 10M:

| Reads | Sensitivity | SNV   | Indel |
|-------|-------------|-------|-------|
| 200K  | 0.272 | 0.381 | 0.141 |
| 500K  | 0.528 | 0.717 | 0.300 |
| 1M    | 0.600 | 0.790 | 0.371 |
| 2M    | 0.613 | 0.800 | 0.388 |
| 5M    | 0.611 | 0.800 | 0.382 |
| 10M   | 0.611 | 0.800 | 0.382 |

At 0.5% VAF the pattern is identical: sensitivity rises steeply to 2M, then is flat.

The hypothesis going in was that more reads would give more alt molecules per target,
reducing coverage-gap failures at the end anchor. The data refutes this.

---

## Diagnostic: Molecule Count vs Read Depth

The first thing to check was whether more reads actually produce more unique molecules.

| Reads | Molecules | Duplex | Duplex % |
|-------|-----------|--------|----------|
| 200K  | 33,993    | 7,166  | 21.1% |
| 500K  | 114,792   | 24,845 | 21.6% |
| 1M    | 226,775   | 47,187 | 20.8% |
| 2M    | 351,341   | 66,552 | 18.9% |
| 5M    | 444,740   | 75,350 | 16.9% |
| 10M   | 466,566   | 76,500 | 16.4% |

**Molecules per read pair:**
- 200K reads: 0.17 mol/read
- 1M reads:   0.23 mol/read
- 2M reads:   0.18 mol/read
- 5M reads:   0.089 mol/read
- 10M reads:  0.047 mol/read

The conversion rate collapses at high depth. From 5M to 10M reads (2× more reads), molecule
count grows by only 5% (444K → 467K). This is UMI space exhaustion.

---

## Root Cause: UMI Space Saturation

Twist UMIs are 5 bases per strand. The canonical pair `min(umi_a + umi_b, umi_b + umi_a)` has
at most `1024 × 1025 / 2 = 524,800` unique values.

Expected occupancy using the coupon-collector model `1 - exp(-n/N)` where N = 524,800:

| Reads | Expected molecules | Expected occupancy |
|-------|-------------------|-------------------|
| 200K  | ~34K              | 6.5% |
| 500K  | ~97K              | 18.5% |
| 1M    | ~179K             | 34.1% |
| 2M    | ~304K             | 58.0% |
| 5M    | ~447K             | 85.3% |
| 10M   | ~497K             | 94.8% |

Observed counts match this curve. At 2M reads, ~58% of UMI space is occupied. At 10M reads,
~95% is occupied. Additional reads at 10M mostly land on already-seen UMI pairs, adding reads
to existing families rather than creating new molecules. Since variant calling is molecule-level
(not read-level), the effective signal does not grow.

**This is the primary cause of the sensitivity plateau.**

---

## Diagnostic: Pathfind Stats at 2M vs 10M Reads (15ng 2% VAF)

To test whether more molecules change graph structure or path enumeration:

| Metric           | 2M reads | 10M reads | Change |
|-----------------|---------|---------|--------|
| Molecules       | 430,786 | 508,985 | +18%   |
| k-mers observed | 882,291 | 1,235,701 | +40%  |
| with_variants   | 277     | 282     | +5     |
| no_paths        | 65      | 62      | -3     |
| ref_only        | 33      | 31      | -2     |
| alt_found       | 277     | 282     | +5     |
| Pathfind time   | 3,387ms | 10,911ms | 3.2×  |
| PASS calls (discovery) | 407 | 541 | +33%  |

The graph has 40% more k-mers at 10M but the pathfind results barely change: only 5 more
targets find alt paths. The extra k-mers are resequencing noise on already-seen molecules, not
new variant-covering paths. More importantly, the 10M run produces 33% more PASS calls in
discovery mode (541 vs 407) but sensitivity in monitoring mode is identical — the extra calls
are background biological FPs and spurious paths, not true positives.

**The second cause of the plateau: the extra molecules at 10M increase graph density, filling
the max_paths=100 budget with more spurious paths and partially offsetting the small gain from
the 5 additional targets.**

---

## Summary of Causes

1. **UMI space saturation** (dominant): at 2M reads, ~58% of all possible 5-base duplex UMI
   pairs are already observed. Additional reads pile onto existing molecules. Molecule count
   grows sub-linearly: 5× reads gives only 18% more molecules. Since calling is molecule-level,
   sensitivity saturates with molecule count, not read count.

2. **Graph density / max_paths budget** (secondary): the extra molecules at 5M and 10M densify
   the de Bruijn graph, introducing more spurious alt paths. The max_paths=100 budget fills with
   noise before the true variant path is found, partially cancelling the small gain from
   additional target coverage.

3. **Structural indel limit** (unchanged): the end-anchor displacement problem for indels is not
   improved by more reads. The end anchor position is fixed by the target design; only the
   relaxed end-anchor fix would improve indel sensitivity regardless of depth.

---

## Implications

- **Deeper sequencing with 5-base UMIs does not improve sensitivity beyond 2M reads.**
  The UMI space is the bottleneck, not coverage.

- **To get more unique molecules, longer UMIs are needed.** 8-base UMIs would give 4^16 / 2 ≈
  2.1B unique pairs — effectively never saturated at any practical depth.

- **The fixes that would actually improve sensitivity at current depth:**
  1. Relaxed end-anchor matching for indels (expected: +20–30 percentage points on indel sensitivity)
  2. Reference-first DFS priority / adaptive max_paths (expected: +5–10 percentage points overall)
  3. Longer UMIs enabling genuine 10M+ read throughput (expected: +20–30 percentage points overall)

- **The read depth sweet spot for 5-base Twist UMIs is 1.5–2M reads.** Below 1M, molecule
  count is too low for reliable low-VAF calling. Above 2M, returns diminish sharply.

---

## Reference

- `hash_partition_umi_grouping.md` — Phase 2 clustering, UMI statistics
- `target_length_investigation.md` — max_paths budget and spurious path analysis
- `benchmarking/results/read_depth/` — full results tables for all six depths
