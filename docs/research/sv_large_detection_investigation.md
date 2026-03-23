# Large SV Detection Investigation (BENCH-VF-009)

## Setup

The sv benchmark suite contains two structural variants per dataset:
- **DUP**: 100bp tandem duplication, junction at pos 500 (target chr1:399-699)
- **INV**: 100bp inversion, junction at pos 900 (target chr1:799-1099)

Both truth VCF entries are large variants. Kam detects only partial alleles at their
junctions. The investigation focuses on why:
- DUP at pos 500 is reliably detected at all VAF levels
- INV at pos 900 is unreliably detected at low VAF and StrandBias-filtered at high VAF

## Evidence

### DUP detection

The DUP creates a 1bp insertion at position 498 in the de Bruijn graph. This short path
is easy to walk and score. The evidence is concentrated: every mutant molecule that
crosses the junction contributes to the same 1-k-mer junction node.

Reported VAF for DUP at various truth VAFs:

| Truth VAF | Reported VAF | Filter |
|-----------|-------------|--------|
| 0.1%      | 1.9%        | PASS   |
| 0.5%      | 3.1%        | PASS   |
| 1.0%      | 5.1%        | PASS   |
| 2.0%      | 9.5%        | PASS   |
| 5.0%      | 16.4%       | PASS   |
| 10.0%     | 30.8%       | PASS   |

The DUP reported VAF is consistently higher than the truth VAF. This is an amplification
artefact: the tandem duplication junction appears in reads from both the left and right
ends of the duplicated region. Molecules spanning the junction from either side contribute
to the same junction insertion, doubling the apparent signal.

### INV detection: two failure modes

**Failure mode 1: LowConfidence at low VAF (0.1%–1.25%)**

The INV path requires walking through ~70 k-mers spanning the 100bp inverted sequence.
The caller uses `min_molecules` across the path as evidence. A single k-mer with low
coverage in the inversion walk bottlenecks the entire path's evidence score.

Reported NALT (molecule evidence) for INV vs truth VAF:

| Truth VAF | NALT | Reported VAF | Filter        |
|-----------|------|-------------|---------------|
| 0.25%     | 0    | —           | not detected  |
| 0.5%      | 2    | 0.19%       | PASS          |
| 0.75%     | 2    | 0.19%       | PASS          |
| 1.0%      | 1    | 0.09%       | LowConfidence |
| 1.25%     | 3    | 0.27%       | PASS          |

NALT is always 1–3 regardless of the truth VAF. At 1% truth VAF, 5000× coverage and ~1000
molecules, we expect ~10 INV molecules. But min_molecules over the 100bp inversion path
drops to 1. This is because:

1. The inversion path starts from the reference sequence, enters the inverted sequence, and
   must walk 100 k-mers through a sequence that only exists in inversion-spanning reads.
2. Each k-mer position in the inversion path is supported only by molecules that fully
   span that exact k-mer position within the inversion.
3. The pathfinder scores by `min_molecules`: 1–3 molecules at the worst-covered k-mer
   position sets the evidence for the entire 100bp path.

The stochastic variation (sometimes PASS, sometimes LowConfidence at the same truth VAF)
is caused by whether the seed produces 2+ molecules vs 1 molecule at the worst k-mer
position in the inversion path.

**Failure mode 2: StrandBias filter at moderate/high VAF (≥2%)**

| Truth VAF | NALT | Reported VAF | Filter     |
|-----------|------|-------------|------------|
| 2.0%      | 12   | 1.1%        | StrandBias |
| 3.0%      | 9    | 0.8%        | StrandBias |
| 5.0%      | 23   | 2.1%        | StrandBias |
| 7.0%      | 27   | 2.4%        | StrandBias |
| 10.0%     | 41   | 3.6%        | StrandBias |

At higher VAF, more molecules support the INV path. But the strand bias filter rejects all
of these calls. The cause is structural: inversion junction reads are inherently strand-biased.

A 100bp inversion creates two junctions:
- Left junction: reads crossing from reference into inverted sequence
- Right junction: reads crossing from inverted sequence back to reference

Junction reads at each end come predominantly from one orientation because:
1. The de Bruijn graph walks the inversion path in one direction
2. Reads that span only the left junction contribute forward-strand support
3. Reads that span only the right junction contribute reverse-strand support
4. Very few reads span the entire 100bp inversion

The result is that the ~100bp inversion path in the de Bruijn graph accumulates evidence
mostly from one strand, failing the Fisher's exact test for strand balance.

At 10% VAF with 41 alt molecules, the INV is consistently found but consistently removed
by StrandBias. This is not a low-signal problem — it is a systematic filter behaviour.

### INS and INVDEL

INS (99bp insertion) and INVDEL (inversion + deletion) show near-zero sensitivity. These
variants require even longer paths in the de Bruijn graph. The min_molecules bottleneck
is more severe (100+ k-mers). Detection only appears at ≥1% VAF and is inconsistent.

## Root Causes

Two independent failure modes prevent reliable large SV detection:

**1. min_molecules path scoring**

For large SVs, the path length is 50–150 k-mers. The minimum over that many positions
is almost always 1–3 molecules, regardless of the true variant allele frequency.
Detection depends on whether the simulation seed places 2+ molecules at every k-mer
position in the path — a lottery at low VAF.

**2. StrandBias filter on inversion junction reads**

Inversion junction reads are structurally strand-biased. At high VAF, where signal is
sufficient, the strand bias filter removes all INV calls. This makes the INV undetectable
at moderate-high VAF even when the de Bruijn graph finds the correct path.

## Implications

1. **DUP detection is an artefact of short path length.** The DUP junction appears as a
   1bp insertion, not a 100bp tandem duplication. Sensitivity is high because the path is
   short and the signal is concentrated. The reported VAF is inflated (~2–3× truth) due to
   double-counting of junction molecules from both ends of the duplication.

2. **INV detection requires architectural changes.** Neither failure mode has a simple fix:
   - min_molecules scoring: changing to mean_molecules or sum would increase false positives
     from noisy single-k-mer paths. A proper fix requires dedicated SV evidence aggregation
     that does not use the SNV/indel scoring model.
   - StrandBias filter: inverting junction reads IS strand-biased. Relaxing the filter for
     structural variant paths would require SV-aware strand bias logic.

3. **Large SV detection needs a dedicated evidence model.** The current approach assigns
   SNV/indel statistical evidence (binomial model on molecule counts at a single position)
   to large SVs detected by path walking. This model is incompatible with large SVs because:
   - The variant allele spans many k-mer positions, each with independent coverage
   - Junction reads are strand-biased by design
   - The path walking produces partial alleles, not full SV sequences

4. **Practical threshold for reliable detection.** Based on these results:
   - SVs detected by path walking can be monitored by position (BENCH-VF-008 fix allows TI
     mode to work), but sensitivity is stochastic below 0.5% VAF.
   - Reliable SV sensitivity (>80%) requires VAF ≥ 0.5% AND the SV to be short enough for
     non-biased junction reads.
   - For clinical use, SV detection below 1% VAF requires a different architecture.

## Follow-up Tasks

These findings suggest two potential follow-ups (not attempted in this investigation):

1. **SV-specific scoring**: Replace min_molecules with a dedicated SV evidence aggregator
   that sums junction molecule evidence across both inversion junctions, rather than taking
   the minimum over the full path. This is a new sub-epic in BENCH-VARFORGE or a new epic.

2. **SV-aware strand bias**: For structural variant paths, compute strand bias on junction-
   spanning molecules only (not the full path). Inversion junction reads from each end have
   opposite strand orientations; combining them makes the overall count balanced.

Both require changes to `kam-pathfind` and `kam-call`. The changes are non-trivial and
are architectural, not bug fixes.
