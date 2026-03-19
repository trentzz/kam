# Titration Benchmark

## Purpose

The titration benchmark tests kam's performance on real clinical-grade reference material with known variant content. It uses the Twist Standard cfDNA Reference Panel v2, which contains 375 QC-passing variants across 84 cancer genes at multiple known VAF levels. This benchmark validates that synthetic performance predicts real-world behaviour and that kam handles the full complexity of real library prep, sequencing noise, and cfDNA fragment characteristics.

## Sample Design

The full titration experiment covers:

- 3 DNA input amounts: 5 ng, 15 ng, 30 ng
- 8 VAF levels: 0%, 0.001%, 0.01%, 0.1%, 0.25%, 0.5%, 1%, 2%
- 3 negative controls (0% VAF, one per DNA input amount)

This gives 24 samples total (21 VAF-positive + 3 negative controls).

### Sample Naming

In all output files, run reports, and papers, use generic sample names of the form `Sample_{ng}ng_VAF_{vaf}` where `{ng}` is the DNA input in nanograms and `{vaf}` is the VAF as a decimal string. Examples: `Sample_5ng_VAF_0.001`, `Sample_30ng_VAF_0.1`, `Sample_15ng_VAF_0`. Do not use actual flowcell IDs, barcode sequences, or internal lab identifiers in any published or committed output.

## Data Locations

These paths are on the local analysis machine (`twz-minipc1`) and are not committed to the repository.

| Resource | Path |
|---|---|
| Raw FASTQ files | `/mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs/` |
| Truth variant set | `/mnt/tzeng-local/tzeng-thesis/titration.probes.QC.pass.tsv` |
| Target sequences | `/mnt/tzeng-local/tzeng-thesis/titration-target-sequences/` |
| Panel BED file | `/mnt/tzeng-local/tzeng-thesis/titration-dedup/bams/panel.bed` |
| Reference genome | `/mnt/tzeng-local/tzeng-thesis/data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` |

FASTQ file SHA-256 checksums are recorded in `docs/benchmarking/data/titration_fastq_checksums.tsv`.

## Truth Set

The truth variant set comes from `titration.probes.QC.pass.tsv`. It contains 375 variants: 205 SNPs and 170 indels. Each row records chromosome, position, REF allele, and ALT allele.

A variant is considered present in a sample if its nominal VAF is greater than 0%. The three 0% VAF samples contain no truth variants; any call in those samples is a false positive.

## Target Sequences and k-mer Parameters

Use the 100 bp flank size from `titration-target-sequences/`. With a default k-mer size of 31, a 100 bp flank on each side of the variant gives adequate anchor k-mers and keeps the index size manageable. Record the exact flank size and k value in every run report.

## Scoring Method

For each sample, match kam calls to truth variants. A call is a true positive if its REF and ALT strings match a truth variant that is expected to be present in that sample at its nominal VAF. A truth variant with no matching call is a false negative. A kam call with no matching truth variant is a false positive.

Compute TP, FP, and FN counts per sample. Then aggregate across samples at the same VAF tier to get per-VAF sensitivity, precision, and F1. Report metrics separately for SNPs and indels.

The scoring script lives at `scripts/benchmarking/score_titration.py`. It takes the kam output VCF, the truth TSV, and the sample VAF as inputs and writes a JSON summary using the same schema as the synthetic benchmark.

### Stratification

Report results stratified by:

- VAF level (primary axis)
- DNA input amount (to quantify input-mass effects on sensitivity)
- Variant class (SNP vs. indel)

Do not stratify by cancer gene or variant position in v0. That analysis is deferred until sensitivity at the variant level is well-characterised.

## Expected Challenges

**Low absolute molecule count at 5 ng input.** At 5 ng of cfDNA and 0.001% VAF, the expected number of variant molecules is very low. Sensitivity in this regime is dominated by library complexity, not by the calling algorithm. Report these data points but note the input-mass constraint explicitly in any table or figure caption.

**UMI collisions.** With a 5 bp UMI (1,024 possible values), collision probability increases with depth. The collision detection logic in `kam-assemble` must handle this correctly. Record the observed collision rate per sample in the QC JSON output.

**Skip base QC.** The 2 bp skip bases between the UMI and the insert are monotemplate. Samples with high skip-base error rates may indicate library prep problems. Record the skip base distribution per sample in the QC JSON output and flag samples that deviate from the expected pattern.
