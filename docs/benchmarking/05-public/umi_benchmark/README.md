# UMI Clustering Benchmark Dataset

## Overview

SRR6794144 is a cfDNA panel sequencing run used as the template dataset in the
Scientific Reports 2025 UMI clustering benchmark. It was used to generate
the simulated spike-in datasets that were benchmarked across eight UMI
clustering tools.

This accession contains real Illumina paired-end reads from a cfDNA panel with
12 bp UMIs embedded in the read name. It provides real UMI sequencing data for
validating kam's molecule assembly before any variant calling.

---

## Provenance

| Property | Detail |
|----------|--------|
| Accession | SRR6794144 |
| BioProject | Not specified in source paper |
| Database | NCBI SRA |
| Access | Open (no application required) |
| Data type | Illumina paired-end FASTQ |
| Sample type | cfDNA from cell-free DNA panel sequencing |
| UMI chemistry | 12 bp simplex UMI; embedded in read name |
| Duplex/simplex | Simplex |
| Platform | Illumina |
| Expected size | ~10 GB compressed (~20 GB uncompressed FASTQ pair) |

---

## UMI Chemistry Notes

The UMI in SRR6794144 is NOT in the Twist 5M2S+T format. The 12 bp UMI is
stored in the read name field of the FASTQ record rather than in the first
bases of R1 or R2.

Before running kam, the UMI must be extracted from the read name and
optionally trimmed to produce reads that match kam's expected input format.
Two approaches:

1. Use `umi_tools extract` or a custom script to move the UMI from the read
   name into the first bases of R1 (synthetic 5M2S+T format with 5 bp prefix).
2. Run kam in a future read-name-UMI mode (not yet implemented).

The spacer/skip bases (2 bp in Twist chemistry) are absent. Set the skip
length to 0 when running kam on this data.

---

## Truth Set

There is no SNV/indel truth set associated with the raw SRR6794144 accession
itself. The truth sets in the 2025 paper belong to simulated datasets derived
from SRR6794144 as a background template, not to the original accession.

This dataset is used exclusively for molecule assembly validation, not for
variant calling accuracy measurement.

---

## Preprocessing Steps

1. Download FASTQ with `fasterq-dump` (see `scripts/download_umi_benchmark.sh`).
2. Inspect the read name format to confirm UMI position (typically
   `@instrument:run:flowcell:lane:tile:x:y:UMI`).
3. Extract the 12 bp UMI from the read name using `umi_tools extract`:

```sh
umi_tools extract \
    --bc-pattern=NNNNNNNNNNNN \
    --stdin SRR6794144_1.fastq.gz \
    --stdout SRR6794144_1.umi.fastq.gz \
    --read2-in SRR6794144_2.fastq.gz \
    --read2-out SRR6794144_2.umi.fastq.gz \
    --log extract.log
```

4. Trim adapters if required (check with FastQC first).
5. Run kam assembly with `--umi-length 12 --skip-length 0`.

---

## Expected Files After Download

```
data/umi_benchmark/
├── SRR6794144_1.fastq.gz    — R1 reads (~5 GB)
└── SRR6794144_2.fastq.gz    — R2 reads (~5 GB)
```

---

## Citation

Calib et al. "Benchmarking UMI clustering tools for accurate detection of
low-frequency variants from deep sequencing." *Scientific Reports* 45, 11773
(2025). https://www.nature.com/articles/s41598-025-33128-x
