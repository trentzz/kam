# Configuration Reference

kam accepts a TOML configuration file via `--config`. All pipeline parameters can be set in
this file. CLI flags always take precedence over config file values, which in turn take
precedence over hard-coded defaults.

## Configuration layering

```
hard-coded defaults
      ↓
  config.toml
      ↓
  CLI flags (always win)
```

A minimal config file requires only the paths that are mandatory for the pipeline to run:

```toml
[input]
r1 = "sample_R1.fq.gz"
r2 = "sample_R2.fq.gz"
targets = "panel_targets.fa"

[output]
output_dir = "results"
```

All other sections and fields are optional. Omitting a field uses the default shown below.

---

## [input]

Input file paths.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `r1` | path | Yes | R1 FASTQ input file (gzip-compressed or plain) |
| `r2` | path | Yes | R2 FASTQ input file (gzip-compressed or plain) |
| `targets` | path | Yes | Target sequences FASTA file. Each entry is one target window. |
| `sv_junctions` | path | No | SV junction sequences FASTA. Required for Inversion and InvDel detection. Not needed for LargeDeletion, TandemDuplication, or NovelInsertion. |
| `target_variants` | path | No | VCF of expected somatic variants for tumour-informed monitoring mode. When set, only calls matching this truth set are marked PASS. |
| `fusion_targets` | path | No | Synthetic fusion target sequences FASTA. Required for Fusion (translocation/gene fusion) detection. Each entry ID must follow the format `{name}__{chromA}:{startA}-{endA}__{chromB}:{startB}-{endB}__fusion`. |

Example:

```toml
[input]
r1 = "sample_R1.fq.gz"
r2 = "sample_R2.fq.gz"
targets = "panel.fa"
sv_junctions = "sv_junctions.fa"
target_variants = "tumour_truth.vcf"
fusion_targets = "fusions.fa"
```

---

## [output]

Output settings.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `output_dir` | path | (required) | Directory for all pipeline outputs. Created if it does not exist. |
| `output_format` | string | `"tsv"` | Output format(s). Comma-separated list of: `tsv`, `csv`, `json`, `vcf`. When multiple formats are given, each is written to a separate file. |
| `qc_output` | path | (none) | Path for a merged QC JSON file. Not currently produced; per-stage QC JSONs are always written to the output directory. |

Example:

```toml
[output]
output_dir = "results"
output_format = "tsv,vcf"
```

---

## [chemistry]

Chemistry and read-structure settings. These define how the UMI, skip, and template regions
are extracted from each read.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `preset` | string | `"twist-umi-duplex"` | Chemistry preset name. Currently informational only; the actual parsing is controlled by `umi_length`, `skip_length`, and `duplex`. |
| `umi_length` | integer | `5` | Length of the random UMI region at the start of each read, in bases. Twist Biosciences: 5 bp. |
| `skip_length` | integer | `2` | Length of the monotemplate spacer/skip region after the UMI, in bases. Twist Biosciences: 2 bp. Set to 0 if no skip region is present. |
| `duplex` | boolean | `true` | Whether the chemistry supports duplex assembly (reads from both strands of the same fragment). Set to `false` for simplex protocols. |
| `min_umi_quality` | integer | `20` | Minimum Phred quality for UMI bases. Read pairs with any UMI base below this quality are dropped before molecule grouping. Set to `0` to disable. |
| `min_template_length` | integer | (none) | Minimum template length in bases. Reads shorter than this after stripping the UMI and skip region are dropped. Comment out to disable. |

The template starts at position `umi_length + skip_length`. For Twist UMI duplex: bases 0–4 are
the UMI, bases 5–6 are the skip, and bases 7+ are the template.

Example for Twist UMI duplex:

```toml
[chemistry]
umi_length = 5
skip_length = 2
duplex = true
min_umi_quality = 20
```

Example for a simplex 12 bp UMI protocol:

```toml
[chemistry]
umi_length = 12
skip_length = 0
duplex = false
min_umi_quality = 20
```

---

## [assembly]

Molecule assembly settings.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `min_family_size` | integer | `1` | Minimum reads per UMI family to be included in consensus building. Setting `1` includes singleton families (single-strand molecules with one read). Setting `2` requires at least two reads per molecule, which eliminates true singletons but reduces sensitivity at low read depth. |

Example:

```toml
[assembly]
min_family_size = 2
```

---

## [indexing]

K-mer indexing settings.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `kmer_size` | integer | `31` | K-mer size. Must be between 1 and 31 (inclusive). Larger k improves specificity (fewer spurious paths) but requires reads that fully span the k-mer window. At k=31, reads shorter than 31 bp produce no k-mers. k=31 is recommended for 100 bp targets at 2M reads. |

Example:

```toml
[indexing]
kmer_size = 31
```

---

## [calling]

Variant calling thresholds. Fields prefixed `sv_` apply to structural variant types
(LargeDeletion, TandemDuplication, Inversion, InvDel, NovelInsertion, Fusion). All other fields
apply to standard variants (SNV, MNV, Insertion, Deletion, Complex).

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `min_confidence` | float | `0.99` | Minimum posterior probability for a PASS call. The posterior is derived from a binomial likelihood ratio between the observed VAF and the background error rate. Calls below this threshold are labelled `LowConfidence`. Lower values increase sensitivity but increase false positives. |
| `strand_bias_threshold` | float | `0.01` | Fisher's exact test p-value cutoff for strand bias. Calls with p below this threshold are labelled `StrandBias`. Genuine somatic variants appear on both strands. Artefacts (oxidative damage, end-repair artefacts) are strand-specific. |
| `min_alt_molecules` | integer | `2` | Minimum molecules supporting the alt allele. A single-molecule call (n_alt=1) is accepted only if duplex-confirmed (see `min_alt_duplex_for_single` in the source). Setting to 3 or higher increases specificity but reduces sensitivity for low-VAF variants. |
| `min_alt_duplex` | integer | `0` | Minimum variant-specific duplex molecules for a PASS call. Counts only duplex molecules whose k-mers overlap the alt path (not the flanking anchor k-mers). Default 0 disables this filter. Set to 1 to require duplex confirmation. Recommended only at duplex fractions ≥15% or depth ≥5M reads. |
| `max_vaf` | float | (none) | Maximum VAF for a PASS call. Calls above this are labelled `HighVaf`. For ctDNA somatic calling, germline heterozygous variants have VAF ≈ 0.5. Setting `0.35` eliminates them while retaining all low-VAF somatic calls. |
| `sv_min_confidence` | float | `0.95` | Minimum posterior probability for SV-type calls. A lower threshold than `min_confidence` is appropriate because a large structural event with 2 supporting molecules is qualitatively stronger evidence than 2 molecules supporting a single-base change. |
| `sv_min_alt_molecules` | integer | `1` | Minimum alt molecule count for SV-type calls. A single-molecule SV call can be meaningful in monitoring mode where the target allele is pre-specified. |
| `sv_strand_bias_threshold` | float | `1.0` | Fisher p-value threshold for strand bias on SV-type variants. Default 1.0 disables the filter. Inversion reads are structurally strand-biased due to directional path walking in the de Bruijn graph; a standard strand bias threshold would flag all inversion calls. |
| `ti_position_tolerance` | integer | `0` | Position tolerance in base pairs for tumour-informed matching. Default 0 requires exact coordinate matching. Use a non-zero value for large SVs where the called position may differ from the truth VCF due to breakpoint ambiguity. |
| `ti_rescue` | boolean | `false` | Enable rescue probing for TI targets that produce no matching call. When set alongside `target_variants`, the k-mer index is queried directly for each undetected TI variant. Results appear with `call_source=RESCUED` or `call_source=NO_EVIDENCE`. |

Example:

```toml
[calling]
min_confidence = 0.99
strand_bias_threshold = 0.01
min_alt_molecules = 2
min_alt_duplex = 0
sv_min_confidence = 0.95
sv_min_alt_molecules = 1
sv_strand_bias_threshold = 1.0
max_vaf = 0.35
target_variants = "tumour_truth.vcf"
ti_position_tolerance = 0
```

---

## [logging]

Logging settings.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `log_dir` | path | (none) | Directory for structured log output. Not currently implemented. |
| `log` | list of strings | (none) | Specific log channels to enable (e.g. `["umi", "family"]`). Not currently implemented. |

---

## [runtime]

Runtime settings.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `threads` | integer | (none) | Number of worker threads. Not currently implemented; the pipeline runs single-threaded. |

---

## CLI flags

All config file fields can be overridden via CLI flags. CLI flags always win.

| CLI flag | Config equivalent |
|----------|------------------|
| `--r1` | `[input] r1` |
| `--r2` | `[input] r2` |
| `--targets` | `[input] targets` |
| `--sv-junctions` | `[input] sv_junctions` |
| `--target-variants` | `[input] target_variants` |
| `--fusion-targets` | `[input] fusion_targets` |
| `--output-dir` | `[output] output_dir` |
| `--output-format` | `[output] output_format` |
| `--chemistry` | `[chemistry] preset` |
| `--min-umi-quality` | `[chemistry] min_umi_quality` |
| `--min-template-length` | `[chemistry] min_template_length` |
| `--min-family-size` | `[assembly] min_family_size` |
| `-k` / `--kmer-size` | `[indexing] kmer_size` |
| `--min-confidence` | `[calling] min_confidence` |
| `--strand-bias-threshold` | `[calling] strand_bias_threshold` |
| `--min-alt-molecules` | `[calling] min_alt_molecules` |
| `--min-alt-duplex` | `[calling] min_alt_duplex` |
| `--max-vaf` | `[calling] max_vaf` |
| `--sv-min-confidence` | `[calling] sv_min_confidence` |
| `--sv-min-alt-molecules` | `[calling] sv_min_alt_molecules` |
| `--sv-strand-bias-threshold` | `[calling] sv_strand_bias_threshold` |
| `--ti-position-tolerance` | `[calling] ti_position_tolerance` |
| `--ti-rescue` | `[calling] ti_rescue` |
| `--threads` | `[runtime] threads` |

---

## Example: Twist UMI duplex, discovery mode

```toml
[input]
r1 = "sample_R1.fq.gz"
r2 = "sample_R2.fq.gz"
targets = "panel.fa"

[output]
output_dir = "results"
output_format = "tsv,vcf"

[chemistry]
umi_length = 5
skip_length = 2
duplex = true
min_umi_quality = 20

[assembly]
min_family_size = 1

[indexing]
kmer_size = 31

[calling]
min_confidence = 0.99
strand_bias_threshold = 0.01
min_alt_molecules = 2
min_alt_duplex = 0
sv_min_confidence = 0.95
sv_min_alt_molecules = 1
sv_strand_bias_threshold = 1.0
```

## Example: Tumour-informed monitoring mode

```toml
[input]
r1 = "ctdna_R1.fq.gz"
r2 = "ctdna_R2.fq.gz"
targets = "panel.fa"
target_variants = "tumour_biopsy.vcf"

[output]
output_dir = "monitoring_results"
output_format = "vcf"

[calling]
max_vaf = 0.35
ti_position_tolerance = 0
```

## Example: SV detection

```toml
[input]
r1 = "sample_R1.fq.gz"
r2 = "sample_R2.fq.gz"
targets = "panel.fa"
sv_junctions = "sv_junctions.fa"
fusion_targets = "fusions.fa"

[output]
output_dir = "sv_results"
output_format = "tsv,vcf"

[calling]
sv_min_confidence = 0.95
sv_min_alt_molecules = 1
sv_strand_bias_threshold = 1.0
```
