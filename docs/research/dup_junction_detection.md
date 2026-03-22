# DUP Junction Detection: Investigation and Fix

## Symptom

The DUP target `chr1:449-649_DUP_100bp` was not detected at any VAF level. The pipeline produced `guided_paths=0` for this target, meaning no alt path was found during graph traversal.

## Diagnostic Data

Graph traversal for the DUP junction failed because the junction k-mer chain broke at J[3]. The first 3 junction k-mers (J[0..2]) were present in the raw graph index (supported by 2–3 molecules at vaf005), but J[3] and beyond were absent. The de Bruijn graph therefore had no edge from J[2] to J[3], and DFS terminated without finding a path to the end anchor.

The root cause of the chain break is geometric: the DUP junction spans alt positions 570–600 in the simulated alt genome. At 0.5% VAF, very few reads cross this region. Of those that do, most are R2 reads (reads from the 3' end of fragments that happen to start close enough to the junction). R2 reads entering the k-mer index as `consensus_rev.sequence` contribute their canonical k-mers correctly, but only 2–3 molecules covered any junction k-mer at vaf005. The raw graph requires both a k-mer AND its suffix to be present for an edge to form; at 2–3 molecules per junction k-mer, coverage is too sparse for a connected chain of 30 k-mers.

## Hypothesis Tested

**Hypothesis 1**: Junction k-mers from `sv_junctions.fa` are absent from the reads entirely.

**Disproved.** The check script (`check_junction_kmers.py`) found hits for k-mers near the junction midpoint, but this script was querying k-mers from the TARGET sequence (not the junction file), so it was actually finding reference k-mers. A deeper analysis confirmed that j-mers from `sv_junctions.fa` ARE present in reads at all VAF levels, but only in small numbers at low VAF.

**Hypothesis 2**: The graph-based DFS cannot find the DUP path because the full junction k-mer chain is disconnected.

**Confirmed.** Even with the full junction chain present in the reads at vaf050, DFS failed at J[3..29] because the sub-budget budget of 300 k-mers is consumed before the DFS reaches the end anchor when starting from a junction k-mer mid-way along the chain.

## Fix: Synthetic Alt Path Construction

The fix bypasses the de Bruijn graph entirely for DUP targets. A synthetic alt path is reconstructed from the target and junction sequences, then scored via the canonical k-mer index.

### Junction File Convention

`sv_junctions.fa` encodes each junction as 60 bases: the last 30 bp of the duplicated region (D) followed by the first 30 bp of D. For the DUP junction:
- Left 30: `target[120..150]` = ref[569..599] = end of D
- Right 30: `target[50..80]` = ref[499..529] = start of D

`synthesize_dup_alt_path` uses `find_subseq` to locate both halves in the target, giving:
- `dup_end = 150` (end of D in target coordinates)
- `dup_start = 50` (start of D = anchor base position)

### Off-by-One Error and Its Fix

The initial implementation constructed the alt sequence as:
```
target[0..dup_end] + D + target[dup_end..]  // WRONG
```

This places the duplication copy at the end of D, producing a junction where D ends and D begins again (D[99] → D[0]). However, the simulated alt genome uses the VCF convention where the copy is inserted right after the anchor base (`dup_start` = ref[499]):

```
alt_genome = ref[0..500] + ref[499..599] + ref[500..]
```

Within the target frame, this is:
```
target[0..51] + target[50..150] + target[51..]  // CORRECT
```

The junction in the real alt genome is D[99] → D[1] (ref[598] → ref[500]), not D[99] → D[0] (ref[598] → ref[499]). These differ by one base: ref[499] = 'C' vs ref[500] = 'A'. As a result, the synthetic path's junction k-mers did not match any k-mer in the canonical index, giving `min_molecules = 0` and a LowConfidence call.

The fix changes the alt sequence construction to:
```rust
let anchor_end = dup_start + 1;
alt_seq = target[0..anchor_end] + target[dup_start..dup_end] + target[anchor_end..]
```

This matches the VCF-derived alt genome and produces correct junction k-mers that are present in the reads.

### Always-Apply Logic

A second issue: the synthetic path was only attempted when `guided.is_empty()`. At vaf010, the graph-based traversal found a chimeric bridge path (a short, low-evidence path through an SNV-like branch), which populated `guided` and suppressed the synthetic path. The fix removes the guard: synthetic DUP paths from `sv_junctions.fa` are always added to the candidate path list, alongside any graph-derived paths.

## Results After Fix

| VAF | TandemDuplication | alt_mol | Filter |
|-----|-------------------|---------|--------|
| 0.5% | Detected | 2 | LowConfidence |
| 1.0% | Detected | 1 | LowConfidence |
| 2.0% | Detected | 4 | PASS |
| 5.0% | Detected | 16 | PASS |

The alt_mol counts reflect the minimum over all 30 junction-crossing k-mers. This is lower than the total number of junction-crossing molecules because reads need to span the full k-mer window, and the last junction k-mer (spanning D[99] → D[1..30]) requires a read to extend 30 bp past the junction.

## What Was Implemented

1. **`synthesize_dup_alt_path`**: Changed alt sequence construction from `target[0..dup_end] + D + target[dup_end..]` to `target[0..dup_start+1] + D + target[dup_start+1..]` (VCF anchor-insertion convention).

2. **Always-apply synthetic path**: Removed the `if guided.is_empty()` guard; synthetic paths from `sv_junctions.fa` are now always evaluated, regardless of whether graph traversal found other paths.

3. **`.min(300).max(1)` → `.clamp(1, 300)`**: Clippy fix in `walk.rs`.

## Open Questions

- The VAF estimates for TandemDuplication are systematically low (1.8% observed at 5% truth). This is because `score_path` uses `min_molecules` over the entire path, and junction k-mers have far fewer supporting molecules than reference k-mers. A better estimator would use the molecule count of the junction-crossing k-mers specifically, not the global minimum. This is a future improvement.

- The LowConfidence calls at vaf005 and vaf010 are expected given the limited junction coverage; they represent correct detections at marginally insufficient evidence. The current `min_alt_molecules` threshold (default 2) would suppress the vaf010 call; this threshold should be reviewed for SV detection.
