# Benchmarking

This directory contains all benchmark material for kam. Each variant class has
its own subfolder with independent data, scripts, results, and documentation.

## Subfolders

### `snvindel/`

Evaluation of SNV and indel detection on the Twist cfDNA Pan-Cancer Reference
Standard v2 titration series (375 truth variants, 24 samples, 5–30 ng input
DNA, VAF 0%–2%). See [`snvindel/README.md`](snvindel/README.md).

### `sv/`

Evaluation of structural variant detection (large deletion, tandem duplication,
inversion) on a 2000 bp synthetic reference using varforge. Covers VAF 0.5%–5%
with Twist duplex UMI chemistry. See [`sv/README.md`](sv/README.md).

## Top-level files

- `.gitignore` — excludes large FASTQ and BAM files from the repository.
