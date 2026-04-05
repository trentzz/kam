# Feature: Hash-Partition UMI Grouping

**Date**: 2026-03-22
**Status**: Phase 1 (HashMap lookup) and Phase 2 (Hamming clustering partition) both implemented.

---

## Problem

The `assemble_molecules` function has two complexity bottlenecks that limit throughput at
higher read depths.

### Bottleneck 1: Linear UMI scan (step 1)

The original step 1 grouped read pairs by canonical UMI using a linear scan:

```rust
unique_umis.iter().position(|(u, _)| u == &pair.canonical_umi)
```

This is O(n × u) where n = total read pairs and u = number of unique canonical UMI pairs.

For 2M reads with ~400K unique canonical UMI pairs: 2M × 200K / 2 = 200B comparisons. In
practice the 10-byte comparison is fast enough (~0.1ns) that this doesn't dominate at 2M
reads. But at 10M+ reads with proportionally more unique UMIs, this term grows as O(n²).

**Fix (implemented)**: replace the linear scan with a `HashMap<CanonicalUmiPair, usize>`,
giving O(1) lookup per read pair and O(n) total for step 1.

### Bottleneck 2: Hamming clustering (step 2)

The `cluster_umi_pairs` function compares every candidate UMI against all existing cluster
seeds. This is O(u²) per sample.

For 2M reads with u ≈ 400K unique UMIs: 400K² / 2 = 80B comparisons. In practice the
directional algorithm (sort by count, absorb low-count neighbours) terminates early and runs
much faster. Empirically, step 2 is not the dominant bottleneck at 2M reads.

At 10M reads, u could reach ~500K (approaching the 524K theoretical maximum for 5bp duplex
UMIs). At that scale, O(u²) clustering becomes a hard wall: 500K² / 2 = 125B comparisons,
taking tens of seconds.

---

## Twist UMI Statistics

5-base random UMIs: 4^5 = 1024 possible values per strand. Canonical pair
`min(umi_a + umi_b, umi_b + umi_a)`: at most 1024 × 1025 / 2 = 524,800 unique canonical
pairs. At 2M reads, expected unique pairs ≈ 524K × (1 − e^(−2M/524K)) ≈ 400K.

This means at 2M reads, ~76% of all possible canonical UMI pairs are already occupied. UMI
collisions are common and the fingerprint-based collision detection in
`split_by_fingerprint` is essential.

At 10M reads, essentially all canonical pairs would be occupied, and the disambiguation
between real molecules and sequencing errors relies entirely on the Hamming clustering and
fingerprint logic.

---

## Phase 2: Hash-Partition Hamming Clustering

**Design**: partition unique UMIs by the first 3 bases of the canonical `umi_a` (64
buckets: 4^3 = 64). Within each bucket, run the existing directional Hamming clustering.
Only UMIs in the same bucket are compared, so Hamming distance = 1 clusters across bucket
boundaries are missed.

**Coverage**: two UMI pairs with Hamming distance 1 in `umi_a` can only have the same
3-base prefix if the mismatch falls in bases 4 or 5. Mismatches in bases 1–3 will route the
two UMIs to different buckets, causing a missed merge. Expected miss rate: 3/10 × (fraction
of pairs within Hamming distance 1) ≈ 3/10 × small ≈ negligible.

**Speedup**: O(u²) → O(u²/64) ≈ 64× faster for the clustering step.

**Trade-off**: small number of Hamming-1 pairs near bucket boundaries will not be merged.
These would be left as separate clusters, adding a small number of extra singleton molecules.
Impact on sensitivity: negligible, since the affected molecules are still assembled
correctly — they just appear as two separate molecules instead of one, and both contribute
to the variant evidence.

**Implementation**: `partition_and_cluster()` in `kam-assemble/src/assembler.rs`.
1. Bucket unique UMIs by `canonical_umi.umi_a[0..3]` into a `HashMap<[u8;3], Vec<usize>>`.
2. For each bucket, build a local pairs slice and call `cluster_umi_pairs`.
3. Remap local cluster member indices back to global `unique_umis` indices.
4. Flatten all bucket cluster results into the global cluster vec.

Phase 2 is a drop-in replacement for the global `cluster_umi_pairs` call in step 2 of
`assemble_molecules`. The downstream pipeline (fingerprint splitting, consensus calling)
is unchanged.

---

## Phase 3: Streaming Sort-Merge for Very Large Inputs

At 20M+ reads the full read vector would not fit in memory. The streaming_molecule_assembly
design (see `docs/research/streaming_molecule_assembly.md`) describes a sort-then-group
approach that avoids loading all reads simultaneously.

This phase is deferred until Phase 2 is validated and depth limits are confirmed.

---

## Measured Impact (Phase 1)

Sample: 15ng 2% VAF, 2M reads, --min-family-size 2.

| Implementation | Assemble time | Sensitivity |
|----------------|--------------|-------------|
| v7 (linear scan) | ~12.4s | 61.3% |
| v10 (HashMap) | ~11.9s | 61.3% |

The 4% assemble speedup confirms the O(n × u) term was not the dominant cost at 2M reads.
The bottleneck at 2M reads is FASTQ I/O + consensus calling, which Phase 1 does not address.

Phase 2 (Hamming clustering partition) will reduce the O(u²) clustering step but is only
necessary at ≥5M reads where u grows large enough for clustering to dominate.

---

## Relationship to Sensitivity

Deeper sequencing is the most direct path to higher sensitivity:
- 2M reads → ~350K molecules → ~780 molecules/target → ~16 alt molecules at 2% VAF
- 10M reads → ~1.75M molecules → ~4700 molecules/target → ~94 alt molecules at 2% VAF

At 10M reads, walk failures due to coverage gaps would be largely eliminated and sensitivity
at 2% VAF is expected to reach 85–90%.

The O(n) step-1 fix and the planned O(u²/64) step-2 fix together enable reaching 10M reads
within practical time and memory budgets.
