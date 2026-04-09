# Investigation: DFS Hangs at High Coverage in ins/sv Discovery

**Date**: 2026-04-09
**Affected jobs**: `ins_ml2_vaf0614_0026` (PIDs 425826, 425866, running 7+ days), plus multiple `ins` and `sv` discovery jobs at high coverage and high VAF
**Finding**: The DFS path walker lacks a total-expansion budget. At high coverage, error k-mers create a near-complete branching factor (~4) and a large dead zone, resulting in `4^132` or more dead-end explorations before termination.

---

## Symptom

`kam run` discovery jobs on `ins` and `sv` variant types get stuck indefinitely. Affected jobs pass assembly (`assembly_qc.json` written) and indexing (`index_qc.json` written) but never produce `calls_discovery.tsv`. The jobs run at 100% CPU for days without progress.

RSS grows slowly over time, from ~200 MB to ~800 MB over one to two days. The hang affects only high-coverage (10,000–17,000x) samples at high VAF (5–14%). Jobs on `indel`, `snv`, and `invdel` variant types complete correctly.

---

## Hypothesis

The DFS in `kam-pathfind/src/walk.rs` has path-length and path-count limits but no limit on total node expansions. At high coverage, sequencing errors create a dense error k-mer graph. Targets with a large `max_path_length` override (e.g. `_maxpath400`) leave a wide gap between the reference path length and the walk budget. The DFS exhausts this gap exploring dead-end branches before it can find or return complete paths.

---

## Root Cause

### Graph construction at high coverage

The `run` pipeline builds the raw graph in `kam/src/commands/run.rs` (lines 162–253) with a two-pass approach. Pass 1 marks molecules as on-target if any k-mer matches the allowlist. Pass 2 inserts all raw k-mers from on-target molecules into `raw_graph_index` with a fixed value of `n_molecules=1` (presence flag only). The graph is then built with `min_molecules=1`, so every k-mer from every on-target molecule enters the graph.

The standalone `pathfind` subcommand uses `min_molecules=2` (pathfind.rs line 103), which filters single-observation k-mers. The `run` pipeline does not.

At Q37 base quality, the per-base error rate is ~0.0002. The probability that a specific single-base substitution k-mer appears in at least one of 16,000 molecules is:

```
1 - (1 - 0.0002)^16000 ≈ 0.96
```

Nearly all three possible substitution k-mers at every reference position are present in the graph. The branching factor approaches 4 at every node. The observed graph for the `ins` target has 76,263 k-mers against 269 reference k-mers — a 280x ratio.

### The dead zone

The `ins_targets.fa` target is 299 bp with a `_maxpath400` suffix, giving 269 reference k-mers and a walk budget of 400.

Error bubbles originate at position `i` and diverge for `k-1 = 30` k-mers, then merge back at position `i+k`. Bubbles starting at positions 0–237 merge back at positions 31–268 — within the target window. These produce complete paths and are valid outputs.

Bubbles starting at positions 238–267 merge back at positions 269–298. These positions are past `end_kmer`. They are dead ends: the DFS enters them and must exhaust their subtree before backtracking. Call this the dead zone.

The DFS is depth-first: after finding the reference path, it backtracks from the deepest stack frame (position 267). The error successors at position 267 enter the dead zone. The remaining walk budget from that point is `400 - 268 = 132`. With a branching factor of 4, the dead-zone search at position 267 alone requires exploring approximately `4^132 ≈ 10^79` paths. There are 30 dead-zone positions (238–267), each with 3 error exits. The DFS cannot reach shallower positions with valid complete paths until every dead-zone branch is exhausted.

### Why other variant types complete

The `snvindel` targets are ~201 bp with `max_path ≈ 221`, giving ~50 k-mers of headroom. At typical snvindel coverage with a branching factor of 1.5–2, the worst-case dead-zone search is approximately `2^51 ≈ 10^15` nodes. This is large but completes within days. The stuck jobs require `4^132`.

For `sv` targets: the INV target uses the default `max_path = 320`, giving only 50 k-mers of headroom — the same regime as snvindel. The DUP target uses `_maxpath400`, giving 130 k-mers of headroom and `4^130` dead-zone cost. Only the DUP target gets stuck.

The relationship is exponential. Each extra k-mer of `max_path_length` beyond the reference path length multiplies the dead-zone search space by the branching factor. Setting `max_path_length = 400` instead of 200 does not double the work — it multiplies it by `4^200 ≈ 10^120`.

The summary across targets:

| Target | Length | maxpath | Ref k-mers | Headroom | Dead-zone cost at B=4 |
|---|---|---|---|---|---|
| `ins_targets.fa` INS | 299 bp | 400 | 269 | 131 | `4^132` |
| `sv_suite_targets.fa` DUP | 300 bp | 400 | 270 | 130 | `4^130` |
| `sv_suite_targets.fa` INV | 300 bp | 320 (default) | 270 | 50 | `4^50` |
| `invdel_targets.fa` INVDEL | 200 bp | 400 | 170 | 230 | `4^230` |
| `snvindel_targets.fa` SNV/INDEL | ~201 bp | ~221 | ~171 | ~50 | `4^50` at typical coverage |

Note on `invdel`: the theoretical dead-zone cost is even larger than `ins`. Invdel completes in practice because the high-VAF, high-coverage combinations that produce B≈4 are less common in the invdel training configs.

### Why RSS grows

Each DFS expansion at `walk.rs` line 156 calls `graph.successors(next).to_vec()`, allocating a fresh `Vec<u64>`. Over billions of backtrack cycles across multiple days, these small allocations (~32 bytes each) fragment the heap. The glibc `malloc` brk watermark does not decrease when small allocations are freed. RSS grows proportionally to total allocation history, not to live set size. The working set of the DFS is bounded (stack + path + visited ≤ ~400 KB), but the allocator's internal free lists touch new pages continuously.

---

## Fix

Four fixes are proposed. Fixes 1 and 2 are implemented. Fixes 3 and 4 are future work.

### Fix 1 (immediate): `max_expansions` budget in `WalkConfig`

Add a `max_expansions: usize` field to `WalkConfig` (`walk.rs` line 50). Increment a counter in the DFS inner loop. When the counter exceeds the budget, return whatever completed paths have been found so far.

A safe budget of 500,000 expansions covers normal cases with 5–15x headroom:
- Reference path: ~269 expansions
- 100 error-bubble paths at ~300 expansions each, plus backtrack steps: ~30,000–100,000 total

At 10^9 expansions per second, a 500,000 budget limits worst-case runtime to ~0.5 ms per target. The stuck jobs currently require `4^132` expansions, which would terminate in microseconds.

Risk: if the budget is hit before the alt path is found, the variant is missed (a false negative). At 5–14% VAF the alt path has strong evidence and the DFS finds it early, well before 500,000 expansions. Risk is low for the specific samples causing the hang.

Add a `walk_budget_exceeded` flag to the discovery QC JSON so early termination is visible in pipeline output.

### Fix 2 (short-term): real molecule counts in the raw graph index

`raw_graph_index` currently stores presence flags only (`n_molecules=1` for all entries, `run.rs` lines 213–214). The `min_molecules` filter in `from_index` therefore cannot prune anything above 1.

Plumb real molecule counts into `raw_graph_index`: count how many distinct consensus molecules contribute each raw k-mer. Then `min_molecules=2` would remove k-mers seen in only a single molecule, matching the behaviour of the standalone `pathfind` subcommand.

This reduces graph density and thus the branching factor. At 16,000x coverage, most error k-mers are seen in 2+ molecules, so this fix alone does not eliminate the problem, but it reduces it substantially for moderate coverage cases.

### Fix 3 (medium-term): iterative-deepening DFS

Replace the single DFS pass with successive passes using increasing path-length budgets: `reference_length + 5`, `reference_length + 10`, `reference_length + 20`, `reference_length + 50`, `reference_length + 100`. Each pass has its own small expansion budget.

This finds short, high-quality paths (single-error bubbles, small indels) in early passes before descending into deep dead-zone branches. Large insertions and duplications are found in later passes, but only after shorter paths are already collected.

### Fix 4 (complementary): per-target wall-clock timeout

Record a start time before each target walk. Check elapsed time every 10,000 expansions. If elapsed exceeds 30 seconds, return current completed paths and set the `walk_budget_exceeded` flag. This acts as a safety net on top of Fix 1.

### Recommended order

1. Fix 1 (`max_expansions`): a five-line change that prevents the current pipeline from stalling. Implement immediately.
2. Fix 2 (real molecule counts): a structural improvement that reduces wasted work at the graph level.
3. Fix 3 (iterative-deepening): addresses the root cause cleanly and replaces `max_expansions` once stable.

---

## Result

Fix 1 and Fix 2 implemented in `feat/ml3-training` (2026-04-09).

- `WalkConfig` now has `max_expansions: usize` (default `0` = unlimited). The DFS counter increments on every stack push and breaks with `budget_exceeded = true` when the limit is reached.
- `walk_paths` and `walk_paths_biased` return `(Vec<GraphPath>, bool)` where the bool is `budget_exceeded`.
- All `kam run` and `kam pathfind` call sites set `max_expansions: 2_000_000`. Sub-walks inside `find_alt_paths_from_reference` use `max_expansions: 300_000`.
- `PathfindQc` gains `n_targets_walk_budget_exceeded: u64`, written to `pathfind_qc.json`.
- `raw_graph_index` now stores real per-molecule counts instead of a sentinel `1`.

Fixes 3 and 4 remain as future work.
