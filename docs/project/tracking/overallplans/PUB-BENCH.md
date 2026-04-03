# PUB-BENCH: Public Dataset Benchmarks

**Status**: todo
**Priority**: high
**Branch**: epic/PUB-BENCH

## Goal

Download, prepare, and run kam on three public datasets that represent
different chemistries and variant types. Produce a cross-dataset results table
suitable for a methods paper. After this epic, kam's performance is validated
on data that did not originate from varforge or the in-house Twist panel.

## Motivation

All current benchmarks use synthetic varforge data. A methods paper requires
at least one external validation dataset. Three datasets are targeted:

1. **UMI Clustering Benchmark** (SRR6794144) — simplex UMI sequencing,
   well-characterised error profile, standard benchmark in UMI literature.
2. **SEQC2 HCC1395** — tumour cell line with a published truth VCF, widely used
   for somatic SNV/indel benchmarking.
3. **COLO829** — melanoma cell line with 68 validated structural variants,
   standard SV benchmark.

Each dataset exercises a different aspect of kam: UMI clustering generalisation,
SNV/indel calling accuracy, and SV detection on real data.

## Design

### UMI Clustering Benchmark (SRR6794144)

Chemistry: simplex, 12 bp UMI, no skip. Requires CHEM-CONFIG (generalised UMI)
to be complete. Run kam in discovery mode; compare molecule count distribution
to published reference.

### SEQC2 HCC1395

Chemistry: standard Illumina paired-end, no UMI (or synthetic UMI injection).
May require a bypass mode (`--no-umi`) or a degenerate single-read-per-molecule
assembly path. Evaluate SNV/indel sensitivity and precision against the SEQC2
truth VCF (Tier 1 variants only).

### COLO829 SV

Chemistry: WGS paired-end, no UMI. Design junction k-mers for the 68 validated
SVs from the published truth set. Run in tumour-informed mode. Evaluate
sensitivity per SV type.

### Cross-dataset table

One row per dataset × mode (discovery / tumour-informed) × variant type.
Columns: sensitivity, precision (where truth is available), N variants evaluated.

## Child tasks

| ID | File | Status |
|----|------|--------|
| PUB-001 | todo/pub_001_dataset_download.md | todo |
| PUB-002 | todo/pub_002_umi_clustering.md | todo |
| PUB-003 | todo/pub_003_seqc2.md | todo |
| PUB-004 | todo/pub_004_colo829.md | todo |
| PUB-005 | todo/pub_005_cross_dataset_table.md | todo |

## Dependencies

- CHEM-CONFIG (required for UMI Clustering Benchmark with 12 bp UMI)
- SV-EXPAND (required for full SV type coverage on COLO829)

## Scope

- `docs/benchmarking/public/` — download scripts, configs, results
- `docs/benchmarking/public/datasets/` — provenance metadata per dataset
- `docs/benchmarking/public/results/` — per-dataset TSV results
- `docs/benchmarking/public/cross_dataset_table.tsv` — combined summary

## Out of scope

- Changes to the kam Rust code
- Downloading full WGS datasets (use subsampled or targeted BAMs)
- Statistical significance testing between datasets
