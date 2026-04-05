# Feature: Pathfind Memory Reduction

## Status

Done. All four fixes implemented and all tests passing.

## Problem

The `pathfind` stage uses 12–14 GB of peak RSS at 1M read pairs on high-input-mass
samples (30ng), while all other stages stay under 1.5 GB. The stage OOMs on most
30ng samples and several 15ng samples at 1M reads. At 200K reads it was already
hitting 4.2 GB on two samples.

### Benchmark evidence (1M reads, 24-sample titration)

| Sample               | Assemble RSS | Index RSS | Pathfind RSS | Notes         |
|---|---|---|---|---|
| 5ng (all VAFs)       | ~1.3 GB      | ~1.2 GB   | 1.3–2.2 GB   | All completed |
| 15ng 0.1% VAF        | 1.1 GB       | 1.0 GB    | 12.6 GB      | OK (slow)     |
| 15ng 2% VAF          | 1.2 GB       | 1.1 GB    | 14.5 GB      | OK (slow)     |
| 15ng 0.25/0.5/1% VAF | ~1.2 GB      | ~1.0 GB   | —            | FAIL (crash)  |
| 30ng (most)          | ~1.2 GB      | ~1.0 GB   | >20 GB       | OOM killed    |
| 30ng 1% VAF          | 1.2 GB       | 1.0 GB    | 7.9 GB       | OK            |

---

## Root Cause Analysis

### 1. BFS queue explosion (dominant cause, estimated 8–14 GB)

`walk_paths` in `walk.rs` uses BFS. Each queue entry is a **full clone** of the
partial path accumulated so far:

```rust
let mut queue: VecDeque<Vec<u64>> = VecDeque::new();
queue.push_back(vec![start_kmer]);

while let Some(partial) = queue.pop_front() {
    for &next in graph.successors(last) {
        let mut new_path = partial.clone();   // full copy every branch
        new_path.push(next);
        queue.push_back(new_path);            // all live branches in queue simultaneously
    }
}
```

**Memory per partial path:** a `Vec<u64>` of length L costs `24 + 8L` bytes.
At depth 35 (halfway through a 70-kmer path) each entry costs ~304 bytes.

**Queue width scales exponentially.** With branching factor B at each node and
path length L, the BFS queue holds up to B^(L/2) entries simultaneously at peak.
For B=3 at 5 positions, B=4 at 2 positions: `3^5 * 4^2 = 3,888` concurrent partial
paths × 304 bytes = ~1.2 MB. Per target, manageable — but the real explosion is
in pathological targets where branches cascade.

**`max_paths=100` does not help.** This limit only counts _completed_ paths (those
that reach `end_kmer`). The queue fills with doomed incomplete paths that will
never complete (dead ends, cycles, too-long branches), and these accumulate without
bound until individually popped and discarded. By the time 100 completed paths are
found, the queue may hold millions of incomplete partial paths.

**Why DFS fixes it:** DFS uses a single path buffer plus a stack of
`(node, sibling_index)` pairs. Memory: O(B × L) = at most a few KB per target,
regardless of branching factor. BFS uses O(B^(L/2) × L) — exponential.

### 2. `max_path_length = 500` is wildly too large

Targets are 100bp. At k=31, the expected path is ~70 k-mers. `max_path_length`
is set to 500 — allowing the walk to explore paths 7× longer than necessary.
This is the main reason the BFS can accumulate so many partial paths: they are
allowed to grow to 500 nodes deep before being pruned. Setting `max_path_length`
to ~150 would cut the maximum BFS depth by 3× and reduce worst-case queue size
by many orders of magnitude.

### 3. Cycle detection is O(path_length)

```rust
if partial.contains(&next) { continue; }  // linear scan every step
```

At path length 70: 70 comparisons per successor per node. This explains the
392-second pathfind wall time for 15ng 0.1% VAF — not a memory issue but a
correctness-performance issue that DFS with a `HashSet` solves simultaneously.

### 4. Predecessors map doubles graph memory

`DeBruijnGraph` stores both `successors` and `predecessors` as
`HashMap<u64, Vec<u64>>`. The predecessors map is never queried during path
walking. It doubles the graph's memory footprint. For a typical 70-node graph
this is ~50 KB — small per target, but unnecessary.

### 5. Index double-materialisation (minor, ~2–4 MB)

`read_bincode` deserialises to `Vec<KmerEntry>`, then `entries_to_hash_index`
builds a `HashMap<u64, MoleculeEvidence>` from it. Both live in memory
simultaneously. At ~26,250 entries × 78 bytes ≈ 2 MB each — negligible compared
to the BFS problem, but wasteful.

### 6. Why 30ng is worse than 5ng

At 5ng, coverage gaps leave some target k-mers absent from the index. The graph
has fewer nodes (perhaps 30–40 of the 70), eliminating branches that run through
uncovered positions.

At 30ng, all 70 target k-mers are observed. The graph is fully connected. Every
possible edge among the 70 target k-mers from suffix/prefix overlap is present.
Targets containing repeated (k-1)-mers (homopolymers, low-complexity regions, or
repetitive sequence) create high-degree nodes that generate many branches. One
pathological target can dominate the entire run's memory footprint.

This is a phase transition: below a certain graph density threshold, BFS stays
manageable; above it, the exponential queue growth overwhelms available memory.
The transition falls between 5ng and 15ng at 1M reads.

---

## Planned Fixes

### Fix 1: Switch BFS to iterative DFS (dominant fix)

Replace the BFS `VecDeque<Vec<u64>>` with an iterative DFS using a single path
buffer and a `(node, sibling_index)` stack. Memory: O(B × L) regardless of
branching factor. A `HashSet<u64>` alongside the path replaces the O(L)
`Vec::contains` cycle detection.

**Files:** `kam-pathfind/src/walk.rs`

### Fix 2: Set max_path_length based on target size (trivial, high impact)

Pass `max_path_length = target_seq.len() / (k - 1) + 50` from the pathfind
command instead of using the hardcoded default of 500. For 100bp targets at k=31,
this gives ~54 — still generous headroom for large indels.

**Files:** `kam/src/commands/pathfind.rs`, `kam-pathfind/src/walk.rs`

### Fix 3: Remove predecessors map from DeBruijnGraph

Remove the `predecessors` field and all code that builds/maintains it. Path
walking only uses `successors`.

**Files:** `kam-pathfind/src/graph.rs`

### Fix 4: Pre-filter low-evidence k-mers before graph construction (strong complement)

In `from_index`, skip k-mers whose molecule count in the index is below a
threshold (e.g. `n_molecules < 2`). This removes PCR error k-mers that create
spurious branches — directly addressing the 30ng branching problem.

**Files:** `kam-pathfind/src/graph.rs`, `kam/src/commands/pathfind.rs`

---

## Implementation Order

1. Fix 3 (remove predecessors) — trivial cleanup, zero risk
2. Fix 2 (max_path_length) — one-liner, immediate win
3. Fix 1 (DFS) — main structural change, eliminates exponential scaling
4. Fix 4 (evidence filtering) — complement to Fix 1

---

## Changes Made

| Fix | Status | Commit | Notes |
|---|---|---|---|
| 1: DFS walk | Done | — | walk.rs: iterative DFS, single path buffer, HashSet cycle detection |
| 2: max_path_length from target | Done | — | pathfind.rs + run.rs: target_seq.len() / (k-1) + 50 per target |
| 3: Remove predecessors | Done | — | graph.rs: predecessors HashMap removed, min_molecules param added |
| 4: Filter low-evidence k-mers | Done | — | min_molecules=2 in from_index call; filters error k-mers before graph build |

---

## Test Plan

- All existing `cargo test` must pass after each fix.
- Re-run titration benchmark at 1M reads after each fix and record peak RSS.
- Target: all samples complete under 4 GB; most under 2 GB.
- Sensitivity must not regress after Fix 4 — verify using comparison plots.
