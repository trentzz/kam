# SEQC2 HCC1395 Somatic Mutation Benchmark

## Overview

The SEQC2 consortium benchmark uses the HCC1395 breast cancer cell line
(tumour) and matched HCC1395BL lymphoblastoid cell line (normal). It is the
most widely used community standard for somatic SNV and indel calling. 63
tumour-normal sequencing pairs were generated across 7 centres; a high-
confidence truth set was constructed by majority-vote consensus of six variant
callers.

---

## Provenance

| Property | Detail |
|----------|--------|
| SRA project | SRP162370 |
| Database | NCBI SRA |
| Access | Open (no application required) |
| Variant types | SNVs, indels |
| UMI chemistry | None |
| Platform | Illumina (primary); also PacBio, Ion Torrent, ONT, 10x |
| Sequencing type | WGS + WES |
| Tumour line | HCC1395 (breast cancer) |
| Normal line | HCC1395BL (matched lymphoblastoid) |

---

## Accessions Downloaded by Script

The download script fetches one representative Illumina WGS tumour-normal
pair to limit disk usage. Full project contains 63 pairs.

| Sample | Accession | Role | Expected size |
|--------|-----------|------|--------------|
| HCC1395 tumour | SRR7890824 | Tumour WGS | ~50 GB |
| HCC1395BL normal | SRR7890827 | Matched normal WGS | ~50 GB |

To download additional replicates, add accessions from SRP162370. Run
`esearch -db sra -query SRP162370 | efetch -format runinfo` to list all runs.

---

## Truth Set

The truth set VCF files are downloaded by the script alongside the FASTQ.

| File | Description | URL |
|------|-------------|-----|
| `HCC1395_truth_SNV.vcf.gz` | High-confidence somatic SNVs | http://bit.ly/2DWuXzP |
| `HCC1395_truth_indel.vcf.gz` | High-confidence somatic indels | http://bit.ly/2PVzzLn |
| NCBI FTP mirror | Both files mirrored on NCBI | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/ |

The truth set uses PASS filter to mark high-confidence calls. Evaluate
kam against PASS variants only.

---

## No-UMI Mode Requirement

SEQC2 data has no UMI. To run kam on this dataset, treat each read pair as a
single molecule (1 read pair = 1 molecule, coverage = 1). This requires a
no-UMI bypass mode in kam that is not yet implemented.

Options when the mode becomes available:

1. `--no-umi` flag: skip UMI extraction and molecule grouping; treat each
   read pair as a distinct molecule.
2. Set a dummy UMI from read position (read name hash); this preserves the
   pipeline structure but degrades to per-read counting.

Until this mode exists, this dataset can be used only to test the variant
calling submodule directly against a reference k-mer index built from
simulated data.

---

## Preprocessing Steps

1. Download FASTQ and truth VCFs (see `scripts/download_seqc2.sh`).
2. Run FastQC on both samples to confirm read quality and adapter presence.
3. Trim adapters with `fastp` or `trimmomatic` if required.
4. Run kam in no-UMI mode (once implemented): `--no-umi`.
5. Evaluate calls against truth VCF using `rtg vcfeval` or `hap.py`.

---

## Expected Files After Download

```
data/seqc2/
├── SRR7890824_1.fastq.gz        — tumour R1 (~25 GB)
├── SRR7890824_2.fastq.gz        — tumour R2 (~25 GB)
├── SRR7890827_1.fastq.gz        — normal R1 (~25 GB)
├── SRR7890827_2.fastq.gz        — normal R2 (~25 GB)
├── HCC1395_truth_SNV.vcf.gz     — SNV truth set
├── HCC1395_truth_SNV.vcf.gz.tbi — index
├── HCC1395_truth_indel.vcf.gz   — indel truth set
└── HCC1395_truth_indel.vcf.gz.tbi
```

Total: ~100 GB for one tumour-normal pair.

---

## Evaluation

Use `rtg vcfeval` (recommended) or `hap.py` with the SEQC2 truth VCF.
Restrict evaluation to PASS variants in the truth set.

```sh
rtg vcfeval \
    --baseline HCC1395_truth_SNV.vcf.gz \
    --calls kam_calls.vcf.gz \
    --template GRCh38.sdf \
    --output rtg_eval_snv/ \
    --sample TUMOR
```

---

## Citation

Xiao et al. "Achieving robust somatic mutation detection with deep learning
models derived from reference data sets of a cancer sample." *Genome Biology*
22, 235 (2021).
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02592-9

SEQC2 consortium site: https://sites.google.com/view/seqc2
NCBI FTP: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/
