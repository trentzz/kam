# COLO829 Somatic SV Truth Set

## Overview

COLO829 is a melanoma cell line with a matched normal COLO829BL. The Cell
Genomics 2022 paper constructed a 68-variant somatic SV truth set by
integrating calls from five sequencing technologies: Illumina HiSeq X Ten,
ONT, PacBio, 10x Genomics linked reads, and Bionano optical mapping. This is
the most established freely available somatic SV benchmark.

Raw Illumina sequencing data are available from ENA with no access
restrictions. The truth VCF and all per-tool call sets are available from
Zenodo.

---

## Provenance

| Property | Detail |
|----------|--------|
| ENA project | PRJEB27698 |
| Access | Open (no application required) |
| Variant types | SVs: 38 DEL, 13 TRA, 7 DUP, 7 INV, 3 INS |
| UMI chemistry | None |
| Platform | Illumina HiSeq X Ten (primary), ONT, PacBio, 10x, Bionano |
| Reference genome | GRCh37 (truth VCF); hg38 liftover available |
| Tumour line | COLO829 (melanoma) |
| Normal line | COLO829BL (matched) |
| Truth set location | Zenodo: https://zenodo.org/records/4716169 |
| GitHub | https://github.com/UMCUGenetics/COLO829_somaticSV |

---

## Truth Set: 68 Somatic SVs

| SV type | Count |
|---------|-------|
| DEL (deletion) | 38 |
| TRA (translocation) | 13 |
| DUP (duplication) | 7 |
| INV (inversion) | 7 |
| INS (insertion) | 3 |
| **Total** | **68** |

All 68 SVs were confirmed by at least two independent technologies. The Zenodo
archive contains:

- `COLO829_truthset_somatic_v4.1.vcf` — GRCh37 truth VCF
- `COLO829_truthset_somatic_hg38.vcf` — hg38 liftover
- Per-tool call sets from GRIDSS, DELLY, Manta, LUMPY, novoBreak, and others
- Validation status annotations

---

## Illumina WGS Accessions

The script downloads the Illumina HiSeq X Ten WGS data. ENA run accessions:

| Sample | ENA run | Role | Expected size |
|--------|---------|------|--------------|
| COLO829 tumour | ERR1341793 | Tumour WGS | ~100 GB |
| COLO829BL normal | ERR1341794 | Matched normal WGS | ~80 GB |

Confirm current accessions at ENA:
https://www.ebi.ac.uk/ena/browser/view/PRJEB27698

---

## No-UMI Mode Requirement

COLO829 data has no UMI. The same no-UMI bypass mode required for SEQC2
applies here. See `seqc2/README.md` for details.

---

## Preprocessing Steps

1. Download Illumina WGS FASTQ and truth VCF (see `scripts/download_colo829.sh`).
2. Run FastQC to confirm read quality.
3. Run kam in no-UMI mode (once implemented).
4. Evaluate SV calls against the truth VCF using Truvari.

---

## Expected Files After Download

```
data/colo829/
├── ERR1341793_1.fastq.gz                    — tumour R1 (~50 GB)
├── ERR1341793_2.fastq.gz                    — tumour R2 (~50 GB)
├── ERR1341794_1.fastq.gz                    — normal R1 (~40 GB)
├── ERR1341794_2.fastq.gz                    — normal R2 (~40 GB)
├── COLO829_truthset_somatic_v4.1.vcf        — GRCh37 truth VCF
└── COLO829_truthset_somatic_hg38.vcf        — hg38 liftover
```

Total: ~180 GB for the Illumina WGS pair.

---

## Evaluation

Use Truvari for SV comparison against the truth VCF.

```sh
truvari bench \
    --base COLO829_truthset_somatic_v4.1.vcf \
    --comp kam_sv_calls.vcf.gz \
    --output truvari_eval/ \
    --refdist 1000 \
    --pctsim 0.7 \
    --pctsize 0.7
```

Note: Truvari defaults are designed for germline SVs. Somatic SV evaluation
may require looser matching parameters, particularly for translocations. Refer
to the COLO829 paper's supplementary for the matching strategy used by the
original authors.

---

## Notes on SV Types and kam Compatibility

| SV type | Notes |
|---------|-------|
| DEL (38) | Directly detectable as k-mer absence/presence signatures |
| TRA (13) | Requires cross-chromosome k-mer linkage; most challenging |
| DUP (7) | Detectable via copy-number elevation in k-mer counts |
| INV (7) | Detectable via strand-reversed k-mer signatures |
| INS (3) | Detectable if insertion sequence is in the index |

Translocations are the most challenging class for alignment-free detection.
The 13 TRA events provide a direct test of whether kam can detect inter-
chromosomal rearrangements.

---

## Citation

van Dijk et al. "A multi-platform reference for somatic structural variation
detection." *Cell Genomics* 2, 100082 (2022).
https://www.sciencedirect.com/science/article/pii/S2666979X22000726

Zenodo truth set archive: https://zenodo.org/records/4716169
GitHub repository: https://github.com/UMCUGenetics/COLO829_somaticSV
