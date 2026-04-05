# Paper Outline: kam — Alignment-Free Variant Detection for Duplex UMI Sequencing

This document defines what the paper should contain and why. It is the reference for
future paper revision cycles. The previous pipeline (HUMID + Jellyfish + km + kmtools)
is included here as context only and should not appear in the paper body.

---

## Context (for author reference only, not in paper)

The thesis (UNSW, Nov 2025) evaluated
HUMID+km on the same Twist STDV2 titration dataset. That pipeline achieved near-zero FPs
through tumour-informed filtering: it only reported km hits matching the exact expected
alt allele for each known target mutation. Background cfDNA variants were invisible by
design. The thesis best SNV sensitivity was 77% and indel sensitivity was 38% at 30ng 2% VAF,
with 375+ truth variants (slightly different truth set). There was no discovery mode.

kam started as the Rust successor to that pipeline, adding molecule-level provenance.
The current paper should not reference HUMID+km as a comparison point — it is background
context for why molecule provenance matters, not a competitor to benchmark against.

---

## Paper Positioning

**One-sentence summary**: kam is an alignment-free, molecule-aware variant detection
pipeline for Twist duplex UMI sequencing that supports both somatic discovery and
tumour-informed monitoring modes.

**Key contribution (what makes this interesting)**:
1. Single-tool replacement for a four-tool chain, with molecule-level evidence preserved end-to-end
2. Tumour-informed monitoring mode (`--target-variants`) achieves precision 1.0 with zero FPs
3. Discovery mode is available for de novo somatic profiling
4. Both modes run in 16–22 seconds per sample on a single core

---

## Sections and Required Content

### Abstract
- One paragraph on the problem (ctDNA, duplex UMI, molecule information loss)
- One paragraph on what kam does (alignment-free, molecule-aware, two modes)
- Numbers: sensitivity 52–61% at 2% VAF, precision 1.0 in monitoring mode, 16–22s runtime
- Mention bug fixes (strand bias, duplex orientation) — these validate the evaluation

### 1. Introduction
- ctDNA clinical motivation (brief)
- Information loss in current pipelines: HUMID discards molecule identity; Jellyfish sees only counts
- Two consequences: PCR duplicate inflation, no duplex/strand distinction at call stage
- kam: single tool, molecule provenance preserved
- Two operating modes: discovery and monitoring
- Key result: monitoring mode → precision 1.0, zero FPs; discovery mode → standard somatic discovery
- Paper structure

### 2. Background
- Duplex UMI sequencing: how it works, error suppression model
- k-mer variant detection: de Bruijn graph, path walking, how it avoids alignment
- Why molecule-level evidence matters (vs read counts)

### 3. Related Work
- HUMID + km pipeline: what it does, what it loses
- Other UMI tools (UMI-tools, fgbio, Sentieon): alignment-based, different design point
- Alignment-based ctDNA callers: Mutect2, TNscope — more sensitive but require alignment
- Do NOT include quantitative comparisons to HUMID+km in this section

### 4. Method
- Molecule assembly: UMI extraction, family grouping, consensus, duplex detection
- k-mer indexing: MoleculeEvidence per k-mer (n_molecules, n_duplex, strand counts)
- de Bruijn graph: construction, DFS path walking, anchor validation
- Variant calling: binomial posterior, strand bias test, quality filters
- Tumour-informed filtering: extract (CHROM, POS, REF, ALT) from called paths, match to target VCF
- Left-normalisation: ensures indel positions match VCF convention
- Recommended configuration: `--max-vaf 0.35 --min-family-size 2 --target-variants` for ctDNA monitoring

### 5. Results
Focus: what does kam achieve on the benchmark dataset?

Required content:
- Dataset: Twist STDV2 24 samples, 375 variants, 3 concentrations, 8 VAF levels, 4 depths
- All reported numbers use: `--max-vaf 0.35 --min-family-size 2 --target-variants truth.vcf`
- Molecule assembly statistics (compression ratio, duplex fraction)
- Sensitivity by VAF and concentration at 2M reads
- Key numbers at 2% VAF: 15ng 61.3%, 30ng 59.2%, 5ng 51.7% — precision 1.0 everywhere
- SNV vs indel sensitivity: SNV 69–80%, indel 31–39%
- Depth scaling (250K to 2M): gains flatten above 1M
- Negative control: 0% VAF → zero calls in monitoring mode
- Runtime: 16–22 s per sample, peak RSS 1.8–2.0 GB

Required figures:
- `sensitivity_vs_vaf.pdf` — main result
- `precision_recall.pdf` — shows precision=1.0
- `confusion_matrix.pdf` — TP/FP/FN at 2% VAF
- `sensitivity_by_type.pdf` — SNV vs indel breakdown
- `compare_sensitivity_vs_vaf.pdf` — depth scaling
- `runtime_per_sample.pdf` — per-sample wall-clock time stacked by stage (assemble/index/pathfind/call) with molecule count per sample as a second panel; shows assembly dominates and runtime scales with molecule yield

### 6. Evaluation
Focus: how did we validate the implementation and what is the operating envelope?

Subsections:
1. **Strand bias overcounting bug** — inflated Fisher test ~70×, how it was found, fix
2. **Duplex orientation comparison bug** — XOR always failed for genuine pairs, fix
3. **Tumour-informed monitoring mode** — the key FP elimination mechanism
   - Discovery mode FPs: 35–72 per sample after statistical filters
   - Why they cannot be removed statistically (overlap with true calls in all evidence dimensions)
   - Monitoring mode: exact allele match to target VCF → 0 FPs
   - Sensitivity unchanged (TPs by definition match the target set)
4. **Sensitivity at sub-0.1% VAF** — why it is near-zero, depth requirement
5. **Missed variant analysis** — FN breakdown by type, indel structural limitation
6. **False positive analysis (discovery mode)** — germline VAF filter, singleton filter, biological floor

### 7. Discussion
- Summary of two-mode design and when to use each
- Tumour-informed filtering as the key FP mechanism: analogy to matched-normal
- SNV vs indel gap and why it persists
- Depth scaling interpretation
- Comparison to alignment-based approaches (qualitative, no quantitative benchmarks)
- Information preservation advantage (molecule evidence for future scoring improvements)
- Reproducibility

### 8. Future Work
- Hash-partition UMI grouping: O(n²) → O(n), enable deeper sequencing
- Head-to-head comparison with Mutect2/TNscope/fgbio+bwa on same read sets
- Duplex-aware posterior scoring (implemented, not yet tuned)
- De novo variant discovery with matched normal annotation
- Nextflow pipeline integration

### Appendix (optional)
- Algorithm pseudocode for path walking
- Detailed parameter table

---

## Numbers Reference (v7 results, 2M reads, tumour-informed mode)

| Sample | VAF | Sensitivity | Precision | SNV sens | Indel sens |
|--------|-----|-------------|-----------|----------|-----------|
| 15ng   | 2%  | 61.3%       | 1.0       | 80.0%    | 38.8%     |
| 30ng   | 2%  | 59.2%       | 1.0       | 77.1%    | 37.6%     |
| 5ng    | 2%  | 51.7%       | 1.0       | 68.8%    | 31.2%     |
| 15ng   | 0.5%| 40.0%       | 1.0       | 52.7%    | 24.7%     |
| 30ng   | 0.5%| 46.1%       | 1.0       | 61.0%    | 28.2%     |
| 5ng    | 0.5%| 16.8%       | 1.0       | 22.4%    | 10.0%     |
| 15ng   | 0.25%| 26.1%      | 1.0       | 36.6%    | 13.5%     |
| 30ng   | 0.25%| 28.5%      | 1.0       | 38.5%    | 16.5%     |
| 5ng    | 0.25%| 6.4%       | 1.0       | 8.3%     | 4.1%      |
| All 0% VAF | — | 0% | — | 0% | 0% (0 FPs) |

Runtime: 16.3–21.7 s per sample. Peak RSS: 1.8–2.0 GB.

---

## What NOT to Include

- Quantitative comparison to HUMID+km (thesis pipeline): it used a different truth set
  version, different mode (monitoring only), and is the predecessor not a competitor
- Detailed HUMID algorithm description: covered in related work only
- Comparison to any alignment-based caller: deferred to future work
- Nextflow workflow details: not yet released
