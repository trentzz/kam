# Synthetic Twist Benchmark

## Purpose

The synthetic benchmark tests whether kam correctly detects known variants at known VAFs in simulated duplex UMI data. Because the ground truth is exact (varforge writes a truth VCF), every false negative and false positive can be attributed to kam's algorithm rather than to uncertainty in the real sample composition. This benchmark isolates algorithmic performance from real-data confounders such as library prep variability, input DNA degradation, and batch effects.

## Data Generation

Synthetic FASTQ pairs are generated with varforge (`/tmp/varforge/target/release/varforge`). The committed config files live in `docs/benchmarking/configs/synthetic/`.

### Key Parameters

| Parameter | Value | Rationale |
|---|---|---|
| UMI structure | 5 bp random + 2 bp skip (`5M2S+T`) | Matches Twist duplex chemistry |
| Fragment model | cfDNA nucleosomal | Realistic fragment length distribution for liquid biopsy |
| Read length | 150 bp | Matches titration sequencing protocol |
| Sequencing depth | 1000x per target | Sufficient molecule count at 0.01% VAF |
| VAF levels | 0%, 0.001%, 0.01%, 0.1%, 0.25%, 0.5%, 1%, 2% | Mirrors titration series |
| RNG seed | Fixed per config file | Required for Nextflow cache compatibility |

### Variant Types

Each config generates a balanced mix of SNPs, insertions (1-10 bp), and deletions (1-10 bp) at each VAF level. Indel size distribution follows the titration truth set proportions (205 SNPs, 170 indels) to keep the synthetic benchmark representative.

### Replicates

Generate at least 5 independent replicates per VAF level by varying the RNG seed. This gives a variance estimate for sensitivity at each VAF. Replicate seeds are listed in the config files; do not use ad-hoc seeds.

## Truth Determination

varforge writes a `truth.vcf` alongside each synthetic FASTQ pair. Every variant in `truth.vcf` is a planted variant with a known allele fraction. This file is the sole source of truth for the synthetic benchmark. Do not modify it after generation.

## Scoring

Matching uses REF + ALT identity only. kam is alignment-free and does not output CHROM:POS coordinates in v0. A call is a true positive if its REF and ALT strings match a variant in `truth.vcf`. A truth variant with no matching call is a false negative. A kam call with no matching truth variant is a false positive.

Compute TP, FP, and FN counts separately for SNPs and indels. Report sensitivity, precision, and F1 per variant class per VAF level.

### Scoring Script

The scoring script lives at `scripts/benchmarking/score_synthetic.py`. It takes the kam output VCF and the varforge `truth.vcf` as inputs and writes a JSON summary. The JSON schema is defined in `docs/research/output_format_specs.md` under the benchmark summary section.

## Advantages Over Real Data

- Unlimited replicates at no sequencing cost.
- Exact ground truth: no ambiguity about whether a variant is truly present.
- Controlled confounders: each parameter can be varied independently (depth, VAF, fragment model, UMI collision rate).
- Fast iteration: run the full benchmark locally without access to sequencing data.

## Known Limitations

- The synthetic fragment model may not capture all sources of real-data noise (e.g., oxidative damage, polymerase errors specific to a sequencing run).
- varforge's UMI collision model uses the analytical approximation described in `docs/research/endpoint_fingerprinting.md`. Real collision rates may differ.
- Results on synthetic data do not transfer directly to real samples. Use the titration benchmark (see `02_titration.md`) to validate that synthetic performance predicts real performance.
