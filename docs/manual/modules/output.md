# kam Pipeline: Outputs

## Overview

Every stage writes at least two output files: the primary data output and a QC JSON. The primary
outputs are the inputs to the next stage. The QC JSON records diagnostic counters for each stage.

---

## Stage outputs

| Stage | Primary output | QC JSON |
|-------|---------------|---------|
| assemble | `molecules.bin` | `assemble_qc.json` |
| index | `index.bin` | `index_qc.json` |
| pathfind | `paths.bin` | `pathfind_qc.json` |
| call | `variants.tsv` (or other format) | `call_qc.json` |

When running `kam run`, all outputs go to `output_dir`. The final variant file is `variants.tsv`
(or `variants.vcf`, etc., depending on `output_format`).

---

## Variant output formats

### TSV (tab-separated values)

Default format. One header line followed by one line per variant call (both PASS and filtered).

**Columns:**

| Column | Type | Meaning |
|--------|------|---------|
| `target_id` | String | Target region identifier (from FASTA header) |
| `variant_type` | String | One of the 11 variant types (see below) |
| `ref_seq` | String | Full reference path sequence |
| `alt_seq` | String | Full alt path sequence |
| `vaf` | Float (6 dp) | VAF point estimate: n_alt / (n_ref + n_alt) |
| `vaf_ci_low` | Float (6 dp) | Lower bound of 95% Beta credible interval |
| `vaf_ci_high` | Float (6 dp) | Upper bound of 95% Beta credible interval |
| `n_molecules_ref` | Integer | Molecules supporting the reference path |
| `n_molecules_alt` | Integer | Molecules supporting the alt path |
| `n_duplex_alt` | Integer | Variant-specific duplex molecules |
| `n_simplex_alt` | Integer | Alt simplex molecules (n_alt - n_duplex_alt) |
| `strand_bias_p` | Float (6 dp) | Fisher's exact test two-tailed p-value |
| `confidence` | Float (6 dp) | Posterior probability the variant is real |
| `filter` | String | Filter label (see below) |

**Note on `ref_seq` and `alt_seq`**: these are the full reconstructed path sequences spanning
the entire target window. The allele representation is not left-normalised in TSV/CSV output.
For VCF output, minimal left-normalised alleles are computed.

### CSV (comma-separated values)

Same columns as TSV with commas as delimiters.

### JSON

A JSON array of objects, one per variant call. Each object has the same fields as the TSV
columns. An empty call list produces `[]`.

Example:
```json
[
  {
    "target_id": "chr17:7674220-7674320",
    "variant_type": "SNV",
    "ref_seq": "ACGT...ACGT",
    "alt_seq": "ACTT...ACGT",
    "vaf": 0.019800,
    "vaf_ci_low": 0.011200,
    "vaf_ci_high": 0.032500,
    "n_molecules_ref": 495,
    "n_molecules_alt": 10,
    "n_duplex_alt": 3,
    "n_simplex_alt": 7,
    "strand_bias_p": 0.842000,
    "confidence": 0.999999,
    "filter": "PASS"
  }
]
```

### VCF (Variant Call Format 4.3)

When the `target_id` encodes a genomic coordinate in `chrN:START-END` format, the VCF writer
produces a properly placed record with minimal left-normalised alleles. For example:

- Target: `chr17:7674220-7674320` (100 bp window, 0-based start)
- Full ref path: 100 bases with `C` at position 50.
- Full alt path: 100 bases with `T` at position 50.
- VCF record: `chr17  7674270  .  C  T  .  PASS  VAF=0.020000;...`

The POS is computed as `start_position + offset_within_path`. For insertions and deletions,
the VCF alleles use minimal left-normalised representation.

If the `target_id` cannot be parsed as a genomic coordinate (e.g. it is a gene name like
`TP53_exon7`), the writer falls back to the alignment-free format: the `target_id` is used as
CHROM and POS is set to 1.

**VCF header meta-information lines:**

```
##fileformat=VCFv4.3
##source=kam-call
##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">
##INFO=<ID=VAF_LO,Number=1,Type=Float,Description="VAF 95% CI lower">
##INFO=<ID=VAF_HI,Number=1,Type=Float,Description="VAF 95% CI upper">
##INFO=<ID=NREF,Number=1,Type=Integer,Description="Molecules supporting REF">
##INFO=<ID=NALT,Number=1,Type=Integer,Description="Molecules supporting ALT">
##INFO=<ID=NDUPALT,Number=1,Type=Integer,Description="Duplex molecules supporting ALT">
##INFO=<ID=NSIMALT,Number=1,Type=Integer,Description="Simplex molecules supporting ALT">
##INFO=<ID=SBP,Number=1,Type=Float,Description="Strand bias Fisher p-value">
##INFO=<ID=CONF,Number=1,Type=Float,Description="Posterior confidence">
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=StrandBias,Description="Strand bias detected">
##FILTER=<ID=LowConfidence,Description="Low posterior confidence">
##FILTER=<ID=LowDuplex,Description="Insufficient duplex support">
##FILTER=<ID=CollisionRisk,Description="UMI collision risk">
##FILTER=<ID=HighVaf,Description="VAF exceeds maximum threshold (likely germline)">
##FILTER=<ID=NotTargeted,Description="Allele not in target-variants set (tumour-informed monitoring mode)">
```

**VCF INFO field summary:**

| Tag | Meaning |
|-----|---------|
| `VAF` | VAF point estimate |
| `VAF_LO` | Lower bound of 95% credible interval |
| `VAF_HI` | Upper bound of 95% credible interval |
| `NREF` | Reference molecule count |
| `NALT` | Alt molecule count |
| `NDUPALT` | Variant-specific duplex molecule count |
| `NSIMALT` | Alt simplex molecule count |
| `SBP` | Strand bias Fisher p-value |
| `CONF` | Posterior confidence |

### Fusion BND output (VCF)

Fusion calls produce VCF BND (breakend) records following VCF 4.3 notation. Each fusion produces
two BND records, one for each partner locus:

```
chr22  23632500  BCR_ABL1_bnd1  C  C[chr9:130854000[  .  PASS  VAF=0.032;...
chr9   130854000 BCR_ABL1_bnd2  A  ]chr22:23632500]A  .  PASS  VAF=0.032;...
```

The partner coordinate is embedded in the ALT field using standard BND bracket notation. Both
records carry the same INFO fields as standard variant calls, plus a `FUSION_PARTNER` INFO tag
pointing to the partner BND record's ID.

---

## Variant types

The `variant_type` field contains one of these 11 values:

| Type | Meaning |
|------|---------|
| `SNV` | Single nucleotide variant |
| `MNV` | Multi-nucleotide variant |
| `Insertion` | Insertion < 50 bp |
| `Deletion` | Deletion < 50 bp |
| `Complex` | Complex rearrangement |
| `LargeDeletion` | Deletion ≥ 50 bp |
| `TandemDuplication` | Tandem duplication ≥ 50 bp |
| `Inversion` | Inversion ≥ 50 bp |
| `Fusion` | Translocation or gene fusion |
| `InvDel` | Inversion with flanking deletion |
| `NovelInsertion` | Novel insertion ≥ 50 bp without sequence similarity to flanking ref |

---

## Filter labels

| Label | Meaning |
|-------|---------|
| `PASS` | All quality filters passed |
| `LowConfidence` | Posterior confidence below threshold, or too few alt molecules |
| `StrandBias` | Alt allele is predominantly on one strand (Fisher p < threshold) |
| `LowDuplex` | Too few variant-specific duplex molecules (when `min_alt_duplex` is set) |
| `HighVaf` | VAF exceeds `max_vaf` threshold (likely germline heterozygous variant) |
| `CollisionRisk` | UMI collision risk too high to trust (not currently assigned) |
| `NotTargeted` | Monitoring mode: allele does not match `target_variants` list |

---

## QC JSON schemas

All QC JSON files are written at the end of each stage. They follow the same top-level structure:

```json
{
  "stage": "stage_name",
  "version": "x.y.z",
  "passed": true,
  ...stage-specific fields...
}
```

### assemble_qc.json

```json
{
  "stage": "molecule_assembly",
  "version": "...",
  "n_read_pairs_processed": 2000000,
  "n_read_pairs_passed": 1983421,
  "n_read_too_short": 142,
  "n_template_too_short": 0,
  "n_low_umi_quality": 16437,
  "n_molecules": 351341,
  "n_duplex": 66552,
  "n_simplex_fwd": 142300,
  "n_simplex_rev": 141200,
  "n_singletons": 1289,
  "n_families_below_min_size": 0,
  "n_umi_collisions_detected": 412,
  "passed": true
}
```

### index_qc.json

```json
{
  "stage": "kmer_indexing",
  "version": "...",
  "kmer_size": 31,
  "n_targets": 375,
  "n_kmers_observed": 882291,
  "n_molecules_indexed": 351341,
  "passed": true
}
```

### pathfind_qc.json

```json
{
  "stage": "pathfind",
  "version": "...",
  "n_targets": 375,
  "no_anchor": 5,
  "no_paths": 65,
  "ref_only": 33,
  "with_variants": 277,
  "passed": true
}
```

### call_qc.json

```json
{
  "stage": "variant_calling",
  "version": "...",
  "n_variants_called": 908,
  "n_pass": 62,
  "n_filtered": 846,
  "passed": true
}
```

---

## Binary intermediate files

The three intermediate binary files (`molecules.bin`, `index.bin`, `paths.bin`) are serialised
using `bincode`. Each file begins with a header record containing:

- `file_type`: a `FileType` enum tag identifying the content.
- `version`: the kam binary version that wrote the file.

This allows format validation and version checking at each stage boundary.

These files are not intended for direct inspection. Use `kam assemble --output molecules.bin`
followed by subsequent stages, or `kam run` to run the full pipeline without manual file
handling.

---

## Per-sample VCF archive

When `run_titration_batch.py` is run with `--save-vcfs DIR`, two VCF files are saved per sample:

| File | Description |
|------|-------------|
| `<name>.monitoring.vcf` | Calls after tumour-informed filter (when `target_variants` was used) |
| `<name>.discovery.vcf` | Calls in discovery mode (a second kam run without `target_variants`) |

The monitoring VCF shows only variants matching the truth set. The discovery VCF shows everything
that passes quality filters, including background biology. Comparing the two reveals the
biological noise floor specific to each sample.
