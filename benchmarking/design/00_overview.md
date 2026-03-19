# Benchmarking Overview

## Hypothesis

Molecule-level k-mer evidence improves sensitivity at low variant allele frequency (VAF) compared to count-based k-mer approaches such as Jellyfish + kmtools. By preserving the identity of the molecule of origin throughout the pipeline, kam can distinguish a true variant supported by multiple independent molecules from a sequencing artefact supported by reads from the same original molecule. This distinction becomes critical below 0.1% VAF, where the absolute variant molecule count is small and artefact suppression determines the detection limit.

## What Is Measured

| Metric | Definition |
|---|---|
| Sensitivity (recall) | TP / (TP + FN) across all truth variants at a given VAF |
| Specificity | TN / (TN + FP) |
| Precision (PPV) | TP / (TP + FP) |
| F1 score | 2 × (precision × recall) / (precision + recall) |
| Detection limit | Lowest VAF level at which sensitivity exceeds 50% |
| Wall time | Elapsed seconds from FASTQ input to VCF output |
| Peak RSS | Maximum resident set size in MB during the run |

All metrics are reported per VAF tier. Detection limit is reported per variant class (SNP, insertion, deletion) where sample size allows.

## What Is Not Measured and Why

**GATK/Mutect2 comparison.** Mutect2 requires read alignment to a reference genome. kam is alignment-free by design. A fair head-to-head comparison would require either aligning kam's output (which destroys the alignment-free property) or running both tools on different representations of the same data (which introduces confounders). Comparative evaluation against alignment-based callers is deferred to a later phase once kam's core claims are validated internally.

**Absolute genomic coordinates.** kam produces variant calls expressed as k-mer evidence rather than CHROM:POS coordinates. Matching to a reference genome for position-based evaluation is not part of the v0 scope. Truth matching uses REF+ALT allele identity only (see individual benchmark documents for the exact matching rule).

**Multi-sample somatic calling.** All benchmarks treat each sample independently. Shared-panel or tumour-normal pair evaluation is future work.

## Baseline and Comparison

kam v0 does not benchmark against an external tool. The primary claim is that sensitivity at VAF ≤ 0.1% is achievable with a molecule-level approach on duplex UMI data. Internal baselines are:

- Per-VAF sensitivity curves showing where the detection limit falls.
- Comparison across DNA input amounts (5 ng, 15 ng, 30 ng) to quantify input-mass effects on sensitivity.
- Comparison between synthetic and real-data benchmarks to assess whether the synthetic model captures real behaviour.

Comparative evaluation against existing tools (HUMID + Jellyfish + kmtools) is planned for a subsequent benchmarking cycle.

## Hardware

All benchmarks run on `twz-minipc1` (local machine). Record the CPU model, core count, RAM, and storage type in each run report alongside the results. This matters because peak RSS and wall time are hardware-dependent.

## Reproducibility

All steps required for reproducibility:

- All RNG seeds are fixed and specified in the varforge config files committed under `docs/benchmarking/configs/`.
- varforge configs and run scripts are committed to the repository. Do not run benchmarks from uncommitted configurations.
- kam is run at a specific git commit hash. Record the hash in every run report.
- Nextflow pipeline scripts pin the exact kam binary version used.
- For the titration benchmark, the raw FASTQ files are not committed (too large), but their SHA-256 checksums are recorded in `docs/benchmarking/data/titration_fastq_checksums.tsv`.
