# SV Benchmark Investigation

## Setup

The SV benchmark suite covers three structural variant types:
- **sv**: DUP (tandem duplication, 100bp) + INV (inversion, 100bp) at positions 500 and 900
- **ins**: INS (insertion, 99bp) at position 399
- **invdel**: INV (100bp inversion) + DEL (20bp) at positions ~300 and ~500

25 VAF levels × 2 replicates × 3 types = 150 datasets. Each run in discovery and
tumour-informed modes.

Note: DEL variants ≥20bp were excluded from the sv suite due to a varforge engine panic
(slice bounds error in engine.rs:393 for any large deletion regardless of position).

## Results Summary

```
type    mode              sens_0010  sens_0025  sens_0050  sens_0100  sens_0200  sens_0500
sv      discovery         0.5        0.5        1.0        0.5        0.5        0.5
sv      tumour_informed   0.0        0.0        0.0        0.0        0.0        0.0
ins     discovery         0.0        0.0        0.0        0.5        0.5        0.0
ins     tumour_informed   0.0        0.0        0.0        0.0        0.0        0.0
invdel  discovery         0.0        0.0        0.0        0.0        0.0        0.5
invdel  tumour_informed   0.0        0.0        0.0        0.0        0.0        0.0
```

## Key Findings

### 1. Tumour-informed sensitivity is always 0% for SVs

**Cause**: Tumour-informed mode filters calls by exact REF/ALT match against the target
VCF. For SVs, the truth VCF contains full alleles (e.g., 100bp insertion sequence). Kam
detects only partial alleles (e.g., 1bp insertion at the first junction). No REF/ALT match
is possible between partial and full alleles.

**Evidence**: At 1% VAF, discovery found a call at pos 498 with ALT=TC (1bp insertion)
while the truth is REF=C, ALT=C+100bp sequence. The TI filter rejects it as NotTargeted.

**Impact**: Tumour-informed mode is not usable for large SVs until either:
a. Kam detects full SV alleles (requires significant pathfinding improvements), or
b. A position-based tumour-informed filter is implemented.

### 2. DUP detection works; INV detection is unreliable

**Evidence**: At 1% VAF, only pos 500 (DUP) is detected (sensitivity 0.5). At 0.5% VAF,
both pos 500 (DUP) and pos 900 (INV) are detected (sensitivity 1.0). The INV detection
fluctuates between VAF levels.

**Comparison (vaf0050 vs vaf0100 for INV at pos 900)**:
- vaf0050: call at pos 901, PASS, VAF=0.19%
- vaf0100: call at pos 901, LowConfidence, VAF=0.09%

Counterintuitively, the detected VAF at 1% is lower than at 0.5%. This suggests the
partial allele for the INV (a small junction fragment in the de Bruijn graph) appears
inconsistently depending on the exact set of reads in each simulation. Path finding for
100bp inversions is inherently fragile: the junction node is shared by few paths, and
coverage noise determines whether it survives path walking.

**Root cause**: 100bp inversion creates a junction k-mer at each end. With 5000x coverage
at 1% VAF = 50 mutant molecules, only ~50 reads cross each junction. The de Bruijn graph
may or may not produce a clear path through the junction depending on noise in that run.

### 3. INS and INVDEL detection is near-zero

For INS (99bp insertion), detection only appears at ≥1% VAF and even then inconsistently.
For INVDEL, detection only appears near 5% VAF.

These larger SVs require the pathfinder to walk longer paths. The `maxpath` override
(maxpath=400 in targets FASTA) allows paths up to 400 k-mers, sufficient in theory for
a 99bp insertion. But in practice, the insertion path competes with reference paths,
and signal from ~25-50 molecules at low VAFs is insufficient to consistently produce PASS calls.

### 4. Junctions removed to prevent pathfind hang

The `--sv-junctions` flag was tested but causes DFS explosion: junction k-mers create
high-connectivity nodes in the de Bruijn graph, making the path walker explore
exponentially many branches without termination.

Running without junctions, kam relies on reference-branching paths only. This is why
SV detection is partial (only small fragments of the SV are reported), not full alleles.

## Implications for Kam Development

1. **Large SV detection needs architectural work.** The current approach of walking
   de Bruijn graph paths was designed for small variants. For large SVs (≥20bp), full
   allele reconstruction requires either:
   - A dedicated SV caller that uses junction evidence differently
   - Indexed junction k-mers without DFS (e.g., direct k-mer graph traversal)
   - Split-read or structural evidence rather than sequence assembly

2. **Tumour-informed mode for SVs requires position-based matching.** The current
   exact REF/ALT filter cannot work for large SVs where alleles are only partially
   reconstructed. A position-based TI filter (match by position ± tolerance) would
   allow monitoring of known SV loci even with partial allele detection.

3. **Small SVs (indels 5-20bp) may be more tractable.** The pathfind approach works
   for SNVs and small indels. Focusing the SV benchmark on deletions/insertions of
   5-15bp would show where the transition from reliable to unreliable detection occurs.
