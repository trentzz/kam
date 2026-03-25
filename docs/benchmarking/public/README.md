# Public Benchmark Datasets

This directory contains documentation and download scripts for the three public
datasets selected to benchmark kam.

All datasets are freely downloadable without data access agreements. Together
they cover kam's three main evaluation goals:

1. Validate UMI molecule assembly on real UMI sequencing data.
2. Validate somatic SNV/indel calling against a community-standard truth set.
3. Validate somatic SV detection against an experimentally confirmed truth set.

---

## Selected Datasets

### 1. UMI Clustering Benchmark — SRR6794144

**Directory:** `umi_benchmark/`
**Download script:** `scripts/download_umi_benchmark.sh`

cfDNA panel FASTQ with 12 bp simplex UMIs. Used as the template dataset in the
Scientific Reports 2025 UMI clustering benchmark. Provides real UMI read
structure for validating molecule assembly independent of variant truth.

Key properties:
- Simplex 12 bp UMI (read-name embedded)
- Deep amplicon sequencing of cfDNA
- Open access via NCBI SRA

Limitations:
- Simplex only; not Twist 5M2S+T duplex chemistry
- No embedded position-based UMI; UMI is in the read name, not the first bases
  of R1/R2 — preprocessing is required to extract it into the expected format
- No associated spike-in truth set for this specific accession (the truth set
  belongs to the simulated derivative datasets from the paper)

**Use in paper:** Confirms that kam's molecule assembly correctly groups reads
by UMI on real cfDNA data, prior to any variant calling.

**Citation:**
Calib et al. "Benchmarking UMI clustering tools for accurate detection of
low-frequency variants from deep sequencing." *Scientific Reports* 45, 11773
(2025). https://www.nature.com/articles/s41598-025-33128-x

---

### 2. SEQC2 HCC1395 Somatic Mutation Benchmark — SRP162370

**Directory:** `seqc2/`
**Download script:** `scripts/download_seqc2.sh`

Multi-centre somatic benchmark using the HCC1395 breast cancer cell line and
matched normal HCC1395BL. The most widely used community standard for somatic
SNV and indel calling. Truth set constructed by majority-vote consensus of six
callers across 63 tumour-normal pairs.

Key properties:
- SNVs and indels; well-characterised high-confidence truth set
- Open access via SRA; truth VCF from NCBI FTP
- No UMI — each read pair is treated as one molecule

Limitations:
- Not UMI data. kam requires a no-UMI mode (or treating each read pair as a
  single molecule) to use this dataset. This mode is not yet implemented.
- Whole-genome/exome, not a targeted panel — k-mer index will be large.
- Benchmarks the variant calling component of kam only; molecule assembly is
  trivial without UMI deduplication.

**Use in paper:** Validates kam's variant calling against a community-standard
truth set without confounding from UMI chemistry differences.

**Citation:**
Xiao et al. "Achieving robust somatic mutation detection with deep learning
models derived from reference data sets of a cancer sample." *Genome Biology*
22, 235 (2021). https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02592-9

---

### 3. COLO829 Somatic SV Truth Set — PRJEB27698

**Directory:** `colo829/`
**Download script:** `scripts/download_colo829.sh`

Multi-platform somatic SV truth set for the COLO829 melanoma cell line and
matched COLO829BL normal. 68 validated somatic SVs confirmed across five
sequencing technologies. The most established freely available somatic SV
benchmark.

Key properties:
- 68 somatic SVs: 38 DEL, 13 TRA, 7 DUP, 7 INV, 3 INS
- Illumina HiSeq X Ten FASTQ available via ENA (no access restrictions)
- Truth VCF available from Zenodo

Limitations:
- No UMI — same no-UMI mode requirement as SEQC2.
- Truth set is for WGS; targeted panel benchmarking requires subsetting to
  panel regions.
- 68 SVs is a small truth set; sensitivity estimates will have wide confidence
  intervals.

**Use in paper:** Validates kam's SV detection against experimentally confirmed
ground truth. The 68 SVs span all major SV classes, testing the full range of
kam's SV calling.

**Citation:**
van Dijk et al. "A multi-platform reference for somatic structural variation
detection." *Cell Genomics* 2, 100082 (2022).
https://www.sciencedirect.com/science/article/pii/S2666979X22000726

---

## Dataset Selection Rationale

Sixteen candidate datasets were evaluated (see
`docs/vision/results/public_datasets.md`). The three above were selected on
the following criteria:

**Open access without application.** Datasets requiring dbGaP, EGA, or
institutional data access agreements (ICGC PCAWG, Hartwig, CODEC, HiDEF-seq)
were excluded to allow immediate, unconditional replication by reviewers.

**FASTQ available.** Alignment-free k-mer analysis requires raw reads.
Datasets available only as aligned BAM (e.g. some ICGC data) were excluded.

**Well-characterised truth set.** The benchmark is only as good as the truth.
SEQC2 (majority-vote consensus, 63 pairs) and COLO829 (five-technology
integration) meet this bar. The UMI dataset has no per-read SNV truth set, but
serves a different purpose (assembly validation).

**Coverage of kam's claims.** The three datasets together test UMI molecule
assembly (dataset 1), somatic SNV/indel calling (dataset 2), and somatic SV
detection (dataset 3).

**Datasets not selected and why:**

| Dataset | Reason not selected |
|---------|-------------------|
| UMI Caller Benchmark (BMC Genomics 2024) | Real cohort data under controlled access; synthetic template SRR10296599 has no SNV truth set |
| cfDNA SNV/Indel Benchmark (Nat Comms 2025) | SRA accession pending at time of writing |
| Duplex Tech Comparison (bioRxiv 2025) | No Twist chemistry; EGA/dbGaP access required for most samples |
| GIAB HG002 small variant | Germline benchmark; not informative for somatic calling |
| GIAB HG002 mosaic | VAF ≥5% only; not representative of ctDNA range |
| SMaHT SNV/indel | Requires WGS-scale compute; no UMI |
| NCI-H2009 SV | ONT only; alignment-free k-mer approach not validated for long reads |
| GIAB HG002 SV | Germline; not somatic |
| SMaHT MIMS SV | WGS-scale; no UMI; infrastructure not yet in place |
| ICGC PCAWG | Controlled access; portal partially decommissioned |
| Hartwig Foundation | Controlled access; application required |
| CODEC | Controlled access (dbGaP) |
| HiDEF-seq | Controlled access (dbGaP); PacBio long-read format |
| Seraseq ctDNA | Physical reference standard; no public FASTQ |
| Horizon Mimix | Physical reference standard; no public FASTQ |

---

## Directory Structure

```
public/
├── README.md                        — this file
├── scripts/
│   ├── download_umi_benchmark.sh    — SRR6794144 via SRA toolkit
│   ├── download_seqc2.sh            — HCC1395 tumour + normal + truth VCF
│   └── download_colo829.sh          — COLO829 Illumina WGS + truth VCF
├── umi_benchmark/
│   └── README.md                    — UMI benchmark dataset notes
├── seqc2/
│   └── README.md                    — SEQC2 dataset notes
└── colo829/
    └── README.md                    — COLO829 dataset notes
```

---

## Prerequisites

All download scripts require either `fasterq-dump` (SRA Toolkit) or `wget`.
Scripts check for these tools and print a clear error if missing.

Install SRA Toolkit: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

Approximate total disk space: ~180 GB (see individual READMEs for breakdown).
