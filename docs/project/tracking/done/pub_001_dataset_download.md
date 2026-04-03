# PUB-001: Dataset Selection and Download Scripts

**Epic**: PUB-BENCH (docs/claudetracking/overallplans/PUB-BENCH.md)
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Write download scripts for all three public datasets and document the data
provenance for each. Done looks like: three shell scripts (one per dataset)
committed to `docs/benchmarking/public/datasets/`, each containing the exact
SRA or HTTP commands to reproduce the download, and a provenance Markdown file
documenting accession numbers, read counts, and intended use.

## Steps

1. Create `docs/benchmarking/public/datasets/` directory.

2. For each dataset, write a download script:

   **UMI Clustering Benchmark (SRR6794144)**
   - Script: `download_umi_clustering.sh`
   - Source: SRA, accession SRR6794144
   - Command: `fastq-dump --split-files --gzip SRR6794144`
   - Document: sample type, UMI length (12 bp), read length, expected read count.

   **SEQC2 HCC1395**
   - Script: `download_seqc2_hcc1395.sh`
   - Source: SEQC2 consortium FTP or SRA (document the exact accession).
   - Choose a subset: 1 FFPE replicate or 1 titration point.
   - Document: truth VCF source (Tier 1 variants only), sample description.

   **COLO829 SV**
   - Script: `download_colo829.sh`
   - Source: SRA or ENA (document accession).
   - Download a subsampled BAM or WGS FASTQ pair.
   - Document: the 68 validated SV truth set source (Lumpy/DELLY consensus or
     published paper supplementary).

3. Write `docs/benchmarking/public/datasets/PROVENANCE.md`:
   - One section per dataset.
   - Fields: accession, source URL, publication reference, download date,
     read count, chemistry, UMI length (if applicable), intended use.

4. Do not run the downloads yet (they may be large and require HPC access).
   The scripts must be correct and runnable, but actual download is deferred
   to the individual per-dataset tasks.

## Notes

- Never commit the raw FASTQ or BAM files to the repo.
- Add `docs/benchmarking/public/datasets/raw/` to `.gitignore`.
- For SEQC2, prefer the IonTorrent or Illumina WES subset (smaller than WGS).
