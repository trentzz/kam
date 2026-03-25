# Public Benchmark Datasets for Alignment-Free Variant Detection

This document lists publicly available datasets that can be used to benchmark
kam, an alignment-free variant detection tool for duplex UMI sequencing. The
focus is datasets with known variant truth sets for SNVs, indels, and
structural variants.

Datasets are grouped by type. Within each group, entries are ordered from most
to least directly applicable to kam's use case (targeted ctDNA panel,
duplex UMI, Twist chemistry).

---

## 1. UMI and ctDNA Sequencing Benchmarks

These datasets most closely match kam's target use case.

### 1.1 UMI Clustering Benchmark (Scientific Reports, 2025)

**Description:** Synthetic and real datasets designed to benchmark eight UMI
clustering tools (AmpUMI, Calib, CD-HIT, Du Novo, Rainbow, Starcode,
UMICollapse, UMI-Tools) for low-frequency variant detection from deep
amplicon sequencing.

| Property | Detail |
|----------|--------|
| Variant types | SNVs, indels |
| VAF range | 0.025% to 10% (simulated); real data at various depths |
| Sequencing depth | 1,000x to 25,000x (simulated) |
| UMI chemistry | Simulated: random 12-bp UMIs in read name. Real: SRR6794144 (12-bp UMI, exact chemistry not specified) |
| Duplex/simplex | Simplex |
| Platform | Illumina |
| FASTQ available | Yes (simulated data generated from SRR6794144 as template) |
| Truth set format | 100 known spike-in variants at specified VAFs |

**Accessions:**
- Reference template: SRR6794144 (NCBI SRA)
- Sample dataset: SRR6794144 (BioProject not specified in paper; original
  cfDNA panel data from NCBI SRA)

**Published comparison:** Calib recommended as best-balanced tool. Full
sensitivity/precision tables across all eight tools and all VAF/depth
combinations published in paper.

**Paper:** [Benchmarking UMI clustering tools for accurate detection of
low-frequency variants from deep sequencing (Scientific Reports,
2025)](https://www.nature.com/articles/s41598-025-33128-x)

---

### 1.2 UMI-Aware Variant Caller Benchmark on ctDNA (BMC Genomics, 2024)

**Description:** Benchmarking Mutect2, VarScan2, shearwater, and DREAMS-vc
using deep targeted sequencing of cfDNA with 9-bp simplex UMIs from 111
colorectal cancer patients. Also uses a synthetic dataset with 303 COSMIC SNV
variants spiked into real cfDNA data from a healthy individual.

| Property | Detail |
|----------|--------|
| Variant types | SNVs (synthetic); somatic SNVs in real samples |
| VAF range | ~0.005% to ~7.5% |
| Sequencing depth | 200x, 450x, 850x (pre-deduplication); mean unique depth ~11,567x (real samples) |
| UMI chemistry | 9-bp simplex UMI; Illumina NovaSeq |
| Duplex/simplex | Simplex |
| Platform | Illumina NovaSeq |
| Panel | Custom 15,396 bp panel covering 12 recurrently mutated CRC genes |
| FASTQ available | Controlled access only (GenomeDK): https://genome.au.dk/library/GDK000009 |
| Truth set format | Synthetic spikes (COSMIC); matched WES for real samples |

**Accessions:**
- Synthetic template: SRR10296599 (cfDNA from healthy individual, NCBI SRA)
- Real cohort: controlled access via https://genome.au.dk/library/GDK000009
  (Frydendahl et al. 2024, same data)
- Metastatic breast cancer samples: SRR15081468, SRR15081470, SRR15081472,
  SRR15081477, SRR15081480, SRR15081482, SRR15081493, SRR15081494 (COMET
  trial, BioProject PRJNA745047)

**Published comparison:** Shearwater with AND consensus best at high precision;
UMIErrorCorrect highest sensitivity at cost of false positives; Mutect2 most
liberal. Full ROC curves in paper.

**Paper:** [Benchmarking UMI-aware and standard variant callers for
low frequency ctDNA variant detection (BMC Genomics,
2024)](https://link.springer.com/article/10.1186/s12864-024-10737-w)

**PMC:** https://pmc.ncbi.nlm.nih.gov/articles/PMC11370058/

---

### 1.3 cfDNA Somatic Variant Calling Benchmark (Nature Communications, 2025)

**Description:** A benchmark defining ~37,000 high-confidence SNVs and ~58,000
indels derived from deep WGS (150x) and WES (2,000x) of patient-matched
cfDNA. Four cancer patients (1 breast, 3 colorectal) with samples at high
(~40%) and ultra-low (~1%) ctDNA tumour fractions. Designed to benchmark nine
somatic callers across ctDNA range.

| Property | Detail |
|----------|--------|
| Variant types | SNVs, indels |
| ctDNA fraction | ~1% to ~40% |
| Sequencing depth | 75x, 150x WGS; 2,000x WES |
| UMI chemistry | No UMI; KAPA Hyper Prep, random 8-mer indices |
| Platform | Illumina NovaSeq 6000 S4 (2×151 bp) |
| FASTQ available | SRA accession pending publication (preprint available) |
| Truth set format | Consensus of 150x + 2,000x patient-matched WGS/WES |

**Accessions:** SRA BioProject number to be confirmed upon publication (listed
as "to be provided during publication process" in preprint).

**Published comparison:** Strelka2 and SMuRF best for SNV discovery at 150x;
VarScan, VarDict, Mutect2 best for clinical genotyping sensitivity; ABEMUS
improved markedly at 2,000x depth.

**Paper:** [Comprehensive benchmarking of methods for mutation calling in
circulating tumor DNA (Nature Communications,
2025)](https://www.nature.com/articles/s41467-025-67842-x)

---

### 1.4 Duplex Sequencing Technology Comparison (bioRxiv, 2025)

**Description:** Benchmarks six duplex sequencing technologies — CODEC,
CompDuplex-seq, HiDEF-seq, NanoSeq (two versions), ppmSeq, and VISTA-seq —
using cord blood DNA, a COLO829/COLO829BL 1:49 tumour-normal mixture, and
tissue homogenates from six human donors. This is the most comprehensive
head-to-head comparison of duplex chemistries.

| Property | Detail |
|----------|--------|
| Variant types | Somatic SNVs |
| VAF range | ~1–2% (COLO829 1:49 mixture) |
| Platform | Illumina NovaSeq (most); PacBio Revio (HiDEF-seq); Ultima UG100 (ppmSeq) |
| UMI/strand chemistry | Varies by method (see table below) |
| Duplex | Yes — all six methods are duplex |
| FASTQ available | Cord blood + neuron samples via EGA: EGAD00001006459; sperm/blood/saliva via dbGaP: phs003716.v1.p1; aristolochic acid dataset: PRJNA1262723 |

**Chemistry detail per method:**

| Method | Platform | Strand barcoding strategy |
|--------|----------|--------------------------|
| NanoSeq-HpyCH4V | Illumina NovaSeq | Reverse-complement UMIs; restriction enzyme fragmentation |
| NanoSeq-MBN | Illumina NovaSeq | Reverse-complement UMIs; mung bean nuclease |
| CompDuplex-seq | Illumina NovaSeq | Tn5-based strand tags |
| VISTA-seq | Illumina NovaSeq X | 16 distinct Tn5 barcode pairs |
| CODEC | Illumina NovaSeq 6000 | Quadruplex physically linked adapters |
| HiDEF-seq | PacBio Revio | Hairpin physically linked adapters |
| ppmSeq | Ultima UG100 | Emulsion-based co-encapsulation |

**Note on relevance to kam:** None of these methods use the Twist 5M2S+T UMI
chemistry. CODEC is closest in concept (full strand linkage) but very
different in implementation. These datasets are useful for understanding
duplex error profiles, not for testing Twist-specific read parsing.

**Paper:** [Benchmarking of duplex sequencing approaches to reveal somatic
mutation landscapes (bioRxiv,
2025)](https://www.biorxiv.org/content/10.64898/2025.12.12.692823v1)

**PMC:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12724167/

---

## 2. Somatic Variant Benchmarks (General)

These datasets use standard (non-UMI) whole-genome or exome sequencing but
have well-characterised truth sets and are widely used as community standards.
They are useful for validating the variant-calling component of kam
independently of the molecule assembly stage.

### 2.1 SEQC2 Somatic Mutation Benchmark (HCC1395)

**Description:** Multi-centre, multi-platform benchmark from the SEQC2
consortium. Tumour (HCC1395 breast cancer) and matched normal (HCC1395BL)
cell lines. 63 tumour-normal pairs from 7 sequencing centres. Truth set of
high-confidence somatic SNVs and indels created by majority-vote consensus of
six callers. The most widely used somatic benchmark in the field.

| Property | Detail |
|----------|--------|
| Variant types | SNVs, indels |
| Sequencing | WGS + WES; Illumina, PacBio, Ion Torrent, ONT, 10x |
| UMI chemistry | No UMI |
| Platform | Multi-platform (Illumina HiSeq/NovaSeq primary) |
| FASTQ available | Yes, open access |
| Truth set format | VCF with PASS filter; annotated with SnpEff/SnpSift |

**Accessions:**
- SRA project: SRP162370
- Normal sample examples: SRR7890827, ERR194147 (ENA)
- Tumour sample examples: SRR7890824, ERR194146 (ENA)
- RNA-seq tumour: BioProject PRJNA635123; normal: PRJNA504037

**Truth set downloads:**
- SNV truth set: http://bit.ly/2DWuXzP
- Indel truth set: http://bit.ly/2PVzzLn
- NCBI FTP: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/

**Site:** https://sites.google.com/view/seqc2

**Published comparison:** Six callers (MuTect2, SomaticSniper, VarDict, MuSE,
Strelka, TNscope) benchmarked; results published in Genome Biology 2021.

**Paper:** [Achieving robust somatic mutation detection with deep learning
models derived from reference data sets of a cancer sample (Genome Biology,
2021)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02592-9)

---

### 2.2 GIAB HG002 Small Variant Benchmark (v4.2.1 / v5.0q)

**Description:** Germline benchmark for a well-characterised Ashkenazi Jewish
individual (HG002/NA24385). The primary benchmark for germline SNVs and
indels. Not somatic, but useful for validating low-level error rates and
confirming correct molecule assembly behaviour on well-understood data.

| Property | Detail |
|----------|--------|
| Variant types | Germline SNVs, indels (v4.2.1); adds SVs in v5.0q |
| UMI chemistry | No UMI |
| Platform | Multi-platform (Illumina 300x, PacBio, ONT, 10x, etc.) |
| FASTQ available | Yes, open access |
| Truth set format | VCF + BED (high-confidence regions) |

**Accessions/downloads:**
- Benchmark files: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/
- Raw Illumina reads: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/

**Site:** https://www.nist.gov/programs-projects/genome-bottle

---

### 2.3 GIAB HG002 Mosaic/Low-Frequency Variant Benchmark (v1.1, 2025)

**Description:** A benchmark of 85 high-confidence subclonal SNVs in the HG002
reference material, with VAF > 5%. The first NIST benchmark designed
specifically for low-frequency (mosaic) variant detection. Directly relevant
to testing low-VAF recall.

| Property | Detail |
|----------|--------|
| Variant types | SNVs (mosaic/subclonal) |
| VAF range | ~5% and above (naturally occurring mosaic variants) |
| UMI chemistry | No UMI |
| Platform | Multi-platform (Illumina primary) |
| FASTQ available | Yes, open access |
| Truth set format | VCF + BED |

**Downloads:**
- FTP: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/mosaic_v1.10/GRCh38/SNV/
- GitHub: https://github.com/usnistgov/giab-HG002-mosaic-benchmark

**Paper:** [Characterisation of subclonal variants in HG002 Genome in a Bottle
reference material as a resource for benchmarking variant callers (Cell
Genomics, 2025)](https://www.cell.com/cell-genomics/fulltext/S2666-979X(25)00360-X)

---

### 2.4 SMaHT SNV/Indel Benchmark — COLO829BLT50 (2025)

**Description:** The SMaHT Network benchmark using the COLO829BL tumour-normal
mixture at 1:49 ratio, producing pseudo-somatic SNVs and indels at VAFs of
~1–2%. Sequenced at 180–500x per replicate across 5 centres on Illumina;
plus PacBio HiFi and ONT long reads. Eleven algorithms evaluated. Focus on
ultra-low allele fraction detection.

| Property | Detail |
|----------|--------|
| Variant types | SNVs, indels |
| VAF range | ~1–2% (cell-line mixture); ~44,000 pseudo-somatic SNVs + ~2,000 indels |
| Sequencing depth | 180–500x Illumina per replicate (5 replicates) |
| UMI chemistry | No UMI |
| Platform | Illumina, PacBio HiFi, ONT |
| FASTQ available | Yes, open access |
| Truth set format | Assembly-derived; orthogonally validated |

**Accessions:** https://data.smaht.org/ (direct download; no SRA submission
required)

**Published comparison:** Call sets were largely discordant across 11 tools;
sensitivity and precision approached ~80% at ≥300x for >2% VAF; long reads
superior for repeat-associated regions.

**Paper:** [Comprehensive benchmarking of somatic single-nucleotide variant
and indel detection at ultra-low allele fractions using short- and long-read
data (2025)](https://pubmed.ncbi.nlm.nih.gov/41278982/)

---

## 3. Structural Variant Benchmarks

### 3.1 COLO829 Somatic SV Truth Set (Cell Genomics, 2022)

**Description:** Multi-platform somatic SV truth set for the COLO829 melanoma
cell line and matched COLO829BL normal. 68 validated somatic SVs derived from
integrating five sequencing technologies. The most established somatic SV
benchmark. Truth set and all raw sequencing data are freely available.

| Property | Detail |
|----------|--------|
| Variant types | SVs: 38 DEL, 13 TRA, 7 DUP, 7 INV, 3 INS |
| Sequencing | Illumina HiSeq X Ten, ONT, PacBio, 10x Genomics, Bionano optical mapping |
| UMI chemistry | No UMI |
| FASTQ available | Yes, open access via ENA |
| Truth set format | VCF (GRCh37 and hg38 liftover) |

**Accessions:**
- ENA project: PRJEB27698 (all raw sequencing; no access restrictions)
- Truth set VCF + per-tool calls: https://zenodo.org/records/4716169
- GitHub: https://github.com/UMCUGenetics/COLO829_somaticSV

**Published comparison:** Multi-tool comparison across all five platforms
published in Cell Genomics 2022.

**Paper:** [A multi-platform reference for somatic structural variation
detection (Cell Genomics,
2022)](https://www.sciencedirect.com/science/article/pii/S2666979X22000726)

---

### 3.2 NCI-H2009 Lung Cancer Somatic SV Dataset (PRJDB10898)

**Description:** Long-read ONT sequencing of NCI-H2009 lung adenocarcinoma and
matched normal BL2009 cell lines. Used in tandem with COLO829 for benchmarking
eight long-read SV callers (Sniffles, cuteSV, Delly, DeBreak, Dysgu, NanoVar,
SVIM, Severus).

| Property | Detail |
|----------|--------|
| Variant types | Somatic SVs (INS, DEL, INV, DUP, TRA) |
| Sequencing | ONT (MinION, GridION, PromethION, R9.4 flow cells) |
| UMI chemistry | No UMI |
| FASTQ available | Yes, open access |
| Truth set | 284 consensus SVs called by majority of tools |

**Accession:** PRJDB10898 (NCBI SRA / DDBJ)

**Paper:** [Benchmarking long-read structural variant calling tools and
combinations for detecting somatic variants in cancer genomes (Scientific
Reports, 2025)](https://www.nature.com/articles/s41598-025-92750-x)

**PMC:** https://pmc.ncbi.nlm.nih.gov/articles/PMC11906795/

---

### 3.3 GIAB HG002 Germline SV Benchmark (v0.6 / v1.0)

**Description:** Germline SV benchmark for HG002, derived from 19 sequence-
resolved calling methods across diverse technologies. The benchmark contains
12,745 isolated, sequence-resolved insertions (7,281) and deletions (5,464),
all ≥50 bp. While germline (not somatic), this is the standard benchmark for
validating SV detection algorithms.

| Property | Detail |
|----------|--------|
| Variant types | INS, DEL (≥50 bp) |
| UMI chemistry | No UMI |
| Platform | Multi-platform |
| FASTQ available | Yes, open access |
| Truth set format | VCF + BED |

**Downloads:**
- v0.6 FTP: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/
- Latest releases: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/

**Paper:** [A robust benchmark for detection of germline large deletions and
insertions (Nature Biotechnology,
2020)](https://www.nature.com/articles/s41587-020-0538-8)

---

### 3.4 SMaHT MIMS SV Benchmark (2025)

**Description:** Synthetic mosaic SV benchmark combining six HapMap cell
lines at controlled mixing ratios (0.5%–83.5%). Contains 34,140 SVs (20,545
INS, 13,595 DEL) ≥50 bp, validated by ddPCR. Total sequencing depth ~2,300x
across Illumina, PacBio HiFi, and ONT. VAFs as low as 0.25%. The most
comprehensive somatic SV benchmark with ultra-low VAF coverage.

| Property | Detail |
|----------|--------|
| Variant types | INS, DEL (≥50 bp); 57.7% in tandem repeats |
| VAF range | 0.25% to germline (≥25%) |
| Sequencing depth | ~1,745x Illumina cumulative; ~380x HiFi; ~180x ONT |
| UMI chemistry | No UMI |
| Platform | Illumina (5 centres), PacBio HiFi (4 replicates), ONT (2 replicates) |
| FASTQ available | Yes, open access |
| Truth set format | Assembly-derived; ddPCR validated; Truvari-harmonised VCF |

**Accessions:** https://data.smaht.org/ (direct download)

**Paper:** [Comprehensive benchmarking of somatic structural variant detection
at ultra-low allele fractions (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458932/)

---

## 4. Controlled-Access and Application-Required Datasets

These datasets require an application or data access agreement but are
publicly available to academic researchers at no cost.

### 4.1 ICGC PCAWG — Pan-Cancer Analysis of Whole Genomes

**Description:** Somatic variant calls (SNVs, indels, SVs, CNVs) from 2,658
whole-genome tumour-normal pairs across 38 cancer types. The largest
catalogue of cancer somatic mutations from WGS. Truth sets are consensus calls
from multiple pipelines. Raw sequencing data available under controlled access.

| Property | Detail |
|----------|--------|
| Variant types | SNVs, indels, SVs, CNVs |
| UMI chemistry | No UMI |
| Platform | Illumina (majority) |
| FASTQ available | Controlled access; application required |
| Truth set format | Consensus VCF per sample |

**Access:**
- Raw data: EGA under EGAC00001000010
- Processed results: https://dcc.icgc.org/pcawg (portal decommissioned; data
  still downloadable via ICGC ARGO)
- ICGC ARGO docs: https://docs.icgc-argo.org/docs/data-access/icgc-25k-data

**Paper:** [Pan-cancer analysis of whole genomes (Nature,
2020)](https://www.nature.com/articles/s41586-020-1969-6)

---

### 4.2 Hartwig Medical Foundation Database

**Description:** WGS tumour-normal pairs from >5,800 metastatic cancer
patients. The largest metastatic WGS database globally. Includes somatic
variant calls (GRIDSS/PURPLE/LINX pipeline), CNVs, SVs, and treatment
response. Not ctDNA; tumour biopsy WGS. Available free of charge to
academic cancer researchers.

| Property | Detail |
|----------|--------|
| Variant types | SNVs, indels, SVs, CNVs |
| UMI chemistry | No UMI |
| Platform | Illumina WGS |
| FASTQ available | Yes (via controlled access request) |
| Truth set | Pipeline consensus; not formally validated against independent truth |

**Access:** Apply via dataaccess@hartwigmedicalfoundation.nl
- Catalogue: https://catalog.hartwigmedicalfoundation.nl
- Documentation: https://hartwigmedical.github.io/documentation/data-access-request-guide.html

---

### 4.3 CODEC Duplex Sequencing Dataset (dbGaP phs003255, Nature Genetics 2023)

**Description:** Duplex sequencing of sperm, blood, and tissue samples using
CODEC (Concatenating Original Duplex for Error Correction). Demonstrated
detection of somatic mutations at genome-wide scale with 1,000-fold accuracy
improvement over standard NGS. Includes liquid biopsy application data.

| Property | Detail |
|----------|--------|
| Variant types | Somatic SNVs, clonal haematopoiesis mutations |
| UMI chemistry | Duplex; physically linked quadruplex adapters; restriction enzyme fragmentation |
| Platform | Illumina NovaSeq 6000 (2×166 or 2×260 cycles) |
| FASTQ available | Controlled access (dbGaP) |
| Truth set | Validated against independent sequencing |

**Access:**
- dbGaP study: phs003255.v1.p1
- Sperm/blood/saliva cohort: dbGaP phs003716.v1.p1

**Paper:** [Single duplex DNA sequencing with CODEC detects mutations with
high sensitivity (Nature Genetics,
2023)](https://www.nature.com/articles/s41588-023-01376-0)

**PMC:** https://pmc.ncbi.nlm.nih.gov/articles/PMC10181940/

---

### 4.4 HiDEF-seq Dataset (dbGaP phs003604, 2023)

**Description:** Single-molecule, long-read duplex sequencing on PacBio that
detects base substitutions on one or both DNA strands. Used to study mutation
spectra in human tissues. Benchmarked against NanoSeq for comparison.

| Property | Detail |
|----------|--------|
| Variant types | Somatic SNVs |
| UMI chemistry | Duplex; hairpin physically linked adapters; long read |
| Platform | PacBio Revio |
| FASTQ/BAM available | Controlled access (dbGaP): raw subreads as PacBio BAM; NanoSeq comparison data as Illumina FASTQ |
| Truth set | Inter-platform validation against NanoSeq |

**Access:** dbGaP study phs003604.v1.p1

---

## 5. Commercial Reference Standards (No Public FASTQ)

These are physical materials with known variant compositions used in wet-lab
validation. Public FASTQ files do not exist but the variant truth sets are
published, so data generated by sequencing them can be benchmarked.

### 5.1 Seraseq ctDNA Mutation Mix v2 (SeraCare)

Purified DNA containing 40 cancer-relevant somatic variants at allele
frequencies from 0.125% to 5%. Validated by digital PCR and NGS. VAF tiers:
0.1%, 0.25%, 2%, and complete mixes. Useful for absolute sensitivity testing
at defined VAFs.

- Product: https://www.seracare.com/Seraseq-ctDNA-Mutation-Mix-v2-AF025-0710-0142/
- Variants: SNVs and short indels; no SVs

### 5.2 Horizon Discovery Mimix Reference Standards

Structural Multiplex (HD753): 18 validated somatic SNVs + copy number
alterations in gDNA. Quantitative Multiplex (HD200): 11 onco-relevant
mutations at 1%–24.5% VAF in FFPE format.

- Product page: https://horizondiscovery.com/en/reference-standards/products/structural-multiplex-reference-standard-gdna
- Variants: SNVs, small indels, CNVs (HD753); SNVs (HD200)

---

## 6. Notes on Alignment-Free Use

None of the above datasets were designed for alignment-free analysis. Key
points for using them with kam:

1. **Read structure compatibility.** Datasets using the Twist 5M2S+T UMI
   chemistry are not yet publicly available at the time of writing (March
   2026). The closest publicly available chemistry is the 9–12 bp simplex UMI
   formats in datasets 1.1 and 1.2.

2. **FASTQ is required.** Datasets where only BAM files are available (some
   ICGC, Hartwig data) cannot be used directly; raw FASTQ is needed for
   alignment-free k-mer analysis.

3. **Simulated data.** For Twist-specific benchmarking before public data
   exists, the varforge simulation framework (already in kam) is the primary
   route. Datasets 1.1 and 1.2 provide validation of the general k-mer and
   molecule assembly approach on real UMI data.

4. **Alignment-based truth sets are still valid.** Truth sets from alignment-
   based pipelines (SEQC2, COLO829, GIAB) remain valid for benchmarking
   kam's variant detection stage. The alignment-free approach changes how
   molecules are assembled, not what constitutes a true positive.

5. **SV detection.** COLO829 (dataset 3.1) provides the best-characterised
   somatic SV truth set with freely downloadable FASTQ. The SMaHT MIMS
   benchmark (dataset 3.4) provides ultra-low VAF SVs but requires WGS-scale
   sequencing depth to use.

---

## Summary Table

| Dataset | Type | VAF range | Duplex/UMI | FASTQ access | Best use for kam |
|---------|------|-----------|------------|--------------|-----------------|
| UMI Clustering Benchmark (1.1) | SNV/indel | 0.025%–10% | Simplex UMI | Open | UMI molecule assembly validation |
| UMI Caller Benchmark (1.2) | SNV | 0.005%–7.5% | Simplex UMI | Controlled/partial open | Caller comparison baseline |
| cfDNA SNV/Indel Benchmark (1.3) | SNV/indel | 1%–40% | No UMI | Pending | WGS caller comparison |
| Duplex Tech Comparison (1.4) | SNV | ~1–2% | Duplex (6 methods) | Partial open | Duplex chemistry understanding |
| SEQC2 HCC1395 (2.1) | SNV/indel | tumour | No UMI | Open | Somatic caller validation |
| GIAB HG002 small variant (2.2) | SNV/indel | germline | No UMI | Open | Assembly correctness check |
| GIAB HG002 mosaic (2.3) | SNV | ≥5% | No UMI | Open | Low-VAF recall benchmark |
| SMaHT SNV/indel (2.4) | SNV/indel | ~1–2% | No UMI | Open | Ultra-low VAF benchmark |
| COLO829 SV truth set (3.1) | SV | somatic | No UMI | Open | SV detection validation |
| NCI-H2009 SV dataset (3.2) | SV | somatic | No UMI | Open | SV tool comparison |
| GIAB HG002 SV (3.3) | INS/DEL | germline | No UMI | Open | SV recall/precision baseline |
| SMaHT MIMS SV (3.4) | INS/DEL | 0.25%–25% | No UMI | Open | Ultra-low VAF SV benchmark |
| ICGC PCAWG (4.1) | SNV/indel/SV | somatic | No UMI | Controlled | Large-scale somatic validation |
| Hartwig Foundation (4.2) | SNV/indel/SV | somatic | No UMI | Controlled | Metastatic cancer validation |
| CODEC dataset (4.3) | SNV | somatic | Duplex | Controlled | Duplex error model comparison |
| HiDEF-seq dataset (4.4) | SNV | somatic | Duplex | Controlled | Long-read duplex comparison |
