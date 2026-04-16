# kam Pipeline: CLI Reference

## Synopsis

```
kam <SUBCOMMAND> [OPTIONS]
```

## Subcommands

| Subcommand | Description |
|-----------|-------------|
| `assemble` | Parse reads and assemble molecules |
| `index` | Build k-mer index from molecules |
| `pathfind` | Walk de Bruijn graph paths |
| `call` | Call variants from scored paths |
| `run` | Run the full pipeline end-to-end |
| `explore` | Interactively inspect bincode data files |
| `models list` | List all built-in ML models |

For most uses, `kam run` is the correct entry point. The individual subcommands are exposed for debugging, benchmarking individual stages, and integration with Nextflow.

---

## kam run

Run the complete pipeline: assemble, index, pathfind, call.

```
kam run --r1 <R1.fastq.gz> --r2 <R2.fastq.gz> --targets <targets.fa> --output-dir <DIR> [OPTIONS]
```

Or with a config file:

```
kam run --config <config.toml> [OPTIONS]
```

### Required arguments

These are required when no `--config` file is provided. When `--config` is given, any of these can be set in the TOML file instead.

---

#### `--r1`

R1 input FASTQ file. Accepts gzip-compressed (`.fq.gz`, `.fastq.gz`) or plain text.

**Example 1: Gzip-compressed input**
```bash
kam run --r1 patient_042_R1.fq.gz --r2 patient_042_R2.fq.gz \
  --targets twist_cfdna_pancancer.fa --output-dir results/
```

**Example 2: Uncompressed input**
```bash
kam run --r1 sample_R1.fastq --r2 sample_R2.fastq \
  --targets panel_v3.fa --output-dir results/
```

---

#### `--r2`

R2 input FASTQ file. Must be paired with `--r1`. Same format requirements as `--r1`.

**Example 1: Standard paired-end run**
```bash
kam run --r1 ctdna_R1.fq.gz --r2 ctdna_R2.fq.gz \
  --targets panel.fa --output-dir results/
```

**Example 2: Absolute paths on a cluster**
```bash
kam run --r1 /data/sequencing/run_2025/sample_01_R1.fq.gz \
  --r2 /data/sequencing/run_2025/sample_01_R2.fq.gz \
  --targets /ref/panels/twist_cfdna_v2.fa --output-dir /scratch/results/sample_01/
```

---

#### `--targets`

Target sequences FASTA file. Each entry defines one target window. Target IDs with the format `chrN:START-END` (0-based, half-open) enable correct VCF coordinate placement.

**Example 1: Panel FASTA**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets twist_cfdna_pancancer.fa --output-dir results/
```

**Example 2: Custom target panel**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets /ref/panels/custom_hotspot_v1.fa --output-dir results/
```

---

#### `--output-dir`

Directory for all pipeline outputs. Created if it does not exist. Stage-specific QC JSON files, intermediate binaries, and final variant calls are all written here.

**Example 1: Relative path**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/patient_042/
```

**Example 2: Absolute path on scratch storage**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir /scratch/kam_results/batch_2025_04/sample_01/
```

---

#### `--config`

Path to a TOML configuration file. When provided, pipeline parameters are loaded from the file first. CLI flags then override individual values. When absent, the four required path flags (`--r1`, `--r2`, `--targets`, `--output-dir`) must appear on the command line.

See the [Configuration Reference](configuration.md) for the full TOML schema.

**Example 1: Config file only**
```bash
kam run --config configs/discovery.toml
```

**Example 2: Config file with CLI overrides**
```bash
kam run --config configs/discovery.toml \
  --r1 new_sample_R1.fq.gz --r2 new_sample_R2.fq.gz \
  --output-dir results/new_sample/
```

**Example 3: Config for tumour-informed mode, override the VCF**
```bash
kam run --config configs/tumour-informed.toml \
  --target-variants patient_042_mutations.vcf
```

---

### Assembly options

#### `--chemistry-override`

Chemistry preset name. Overrides the config file value. Currently only `twist-umi-duplex` is supported, which sets UMI length = 5 bp, skip length = 2 bp, and template start = 7 bp.

**Example 1: Explicit Twist chemistry**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --chemistry-override twist-umi-duplex
```

**Example 2: Override config file chemistry**
```bash
kam run --config configs/discovery.toml \
  --chemistry-override twist-umi-duplex
```

---

#### `--min-umi-quality-override`

Minimum Phred quality for UMI bases. Read pairs with any UMI base below this quality are dropped before grouping. Default: 20. Set to 0 to disable quality filtering. Overrides the config file value.

**Example 1: Increase stringency to Q30**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-umi-quality-override 30
```

**Example 2: Disable UMI quality filtering**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-umi-quality-override 0
```

**Example 3: Low-quality library, relax threshold**
```bash
kam run --config configs/discovery.toml \
  --min-umi-quality-override 10
```

---

#### `--min-family-size-override`

Minimum reads per UMI family. Families with fewer reads are discarded. Default: 1. Setting to 2 eliminates true singletons but reduces sensitivity at low depth. Overrides the config file value.

**Example 1: Require at least 2 reads per family**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-family-size-override 2
```

**Example 2: Include singletons (default)**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-family-size-override 1
```

**Example 3: High-depth library, stricter families**
```bash
kam run --config configs/high-specificity.toml \
  --min-family-size-override 3
```

---

#### `--min-template-length`

Minimum template length in bases. Templates shorter than this on either R1 or R2 are dropped. Not set by default. Useful for filtering very short inserts that produce low-quality consensus.

**Example 1: Drop templates shorter than 30 bp**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-template-length 30
```

**Example 2: Filter adapter dimers**
```bash
kam run --config configs/discovery.toml \
  --min-template-length 50
```

---

### Indexing options

#### `-k` / `--kmer-size-override`

K-mer length for indexing. Must be 1 to 31. Larger k improves specificity (fewer spurious paths) but requires reads that fully span the k-mer window. At k=31, reads shorter than 31 bp produce no k-mers. Default: 31. Overrides the config file value.

**Example 1: Default k=31 for standard panels**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  -k 31
```

**Example 2: Shorter k for short target windows**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets short_amplicons.fa --output-dir results/ \
  -k 21
```

**Example 3: Override config for benchmarking**
```bash
kam run --config configs/discovery.toml --kmer-size-override 25
```

---

### Calling options

#### `--min-confidence`

Minimum posterior probability for a PASS call. Derived from a binomial likelihood ratio between the observed VAF and the background error rate. Calls below this threshold are labelled `LowConfidence`. Default: 0.99. Lower values increase sensitivity at the cost of more false positives.

**Example 1: Default high-confidence calling**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-confidence 0.99
```

**Example 2: Relaxed for discovery in low-depth samples**
```bash
kam run --r1 low_depth_R1.fq.gz --r2 low_depth_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-confidence 0.90
```

**Example 3: Ultra-stringent for clinical reporting**
```bash
kam run --config configs/tumour-informed.toml \
  --min-confidence 0.999
```

---

#### `--strand-bias-threshold`

Fisher's exact test p-value cutoff for strand bias. Calls with p below this threshold are labelled `StrandBias`. Genuine somatic variants appear on both strands. Artefacts (oxidative damage, end-repair artefacts) tend to be strand-specific. Default: 0.01.

**Example 1: Default strand bias filter**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --strand-bias-threshold 0.01
```

**Example 2: Stricter filter for FFPE samples**
```bash
kam run --r1 ffpe_R1.fq.gz --r2 ffpe_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --strand-bias-threshold 0.05
```

**Example 3: Relaxed filter for low-depth monitoring**
```bash
kam run --config configs/tumour-informed.toml \
  --strand-bias-threshold 0.001
```

---

#### `--min-alt-molecules`

Minimum molecules supporting the alt allele for a PASS call. A single-molecule call (n_alt=1) is accepted only if duplex-confirmed. Default: 2. Setting to 3 or higher increases specificity but reduces sensitivity for low-VAF variants.

**Example 1: Default (2 molecules)**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-alt-molecules 2
```

**Example 2: High specificity for clinical reporting**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-alt-molecules 3
```

**Example 3: Monitoring mode, allow single-molecule calls**
```bash
kam run --r1 ctdna_R1.fq.gz --r2 ctdna_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --target-variants tumour_mutations.vcf \
  --min-alt-molecules 1
```

---

#### `--min-alt-duplex`

Minimum variant-specific duplex molecules. Counts only duplex molecules whose k-mers overlap the variant site (alt-path-specific k-mers), not flanking anchor k-mers. Default: 0 (disabled). Set to 1 to require at least one duplex confirmation at the variant site. Recommended only at duplex fractions above 15% or depth above 5M reads.

**Example 1: Require duplex confirmation**
```bash
kam run --r1 high_depth_R1.fq.gz --r2 high_depth_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-alt-duplex 1
```

**Example 2: Disabled (default)**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-alt-duplex 0
```

**Example 3: Stringent duplex for high-depth clinical run**
```bash
kam run --r1 clinical_R1.fq.gz --r2 clinical_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-alt-duplex 2 --output-format-override tsv,vcf
```

---

#### `--max-vaf`

Maximum VAF for a PASS call. Calls with VAF above this are labelled `HighVaf`. Not set by default. For ctDNA somatic calling, germline heterozygous variants have VAF around 0.5. Setting `--max-vaf 0.35` eliminates them while retaining all low-VAF somatic calls.

**Example 1: Exclude germline variants**
```bash
kam run --r1 ctdna_R1.fq.gz --r2 ctdna_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --max-vaf 0.35
```

**Example 2: Allow higher VAF for post-treatment monitoring**
```bash
kam run --r1 post_chemo_R1.fq.gz --r2 post_chemo_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --target-variants tumour_mutations.vcf \
  --max-vaf 0.45
```

**Example 3: Very low VAF cap for ultra-pure somatic detection**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --max-vaf 0.10
```

---

#### `--target-variants`

VCF file of expected somatic variants for tumour-informed monitoring mode. When set, only calls whose `(CHROM, POS, REF, ALT)` matches an entry in this VCF are marked PASS. All other quality-passing calls are labelled `NotTargeted`. Use this when the somatic variant panel is known from a prior tissue biopsy. Produces near-zero false positives.

**Example 1: Standard tumour-informed monitoring**
```bash
kam run --r1 plasma_tp2_R1.fq.gz --r2 plasma_tp2_R2.fq.gz \
  --targets twist_cfdna_pancancer.fa --output-dir results/tp2/ \
  --target-variants tumour_biopsy_snvs.vcf \
  --max-vaf 0.35
```

**Example 2: Combined with TI rescue for complete tracking**
```bash
kam run --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --target-variants patient_042_mutations.vcf \
  --ti-rescue \
  --output-format-override tsv,vcf
```

**Example 3: Monitoring with position tolerance for SVs**
```bash
kam run --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --target-variants patient_sv_mutations.vcf \
  --ti-position-tolerance-override 10
```

---

#### `--ti-position-tolerance-override`

Position tolerance in base pairs for tumour-informed matching. Default: 0 (exact coordinate matching). When greater than 0, a call also passes the TI filter if its position is within this many bp of any target variant position, regardless of REF/ALT. Useful for large SVs with breakpoint ambiguity where the called position may differ from the truth VCF.

**Example 1: Exact matching (default)**
```bash
kam run --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --target-variants tumour_mutations.vcf \
  --ti-position-tolerance-override 0
```

**Example 2: Allow 1 bp tolerance for indel left-alignment differences**
```bash
kam run --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --target-variants tumour_mutations.vcf \
  --ti-position-tolerance-override 1
```

**Example 3: Large tolerance for SV breakpoints**
```bash
kam run --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --target-variants patient_sv_mutations.vcf \
  --sv-junctions sv_junctions.fa \
  --ti-position-tolerance-override 10
```

---

#### `--ti-rescue`

Enable rescue probing for TI targets that produce no matching call. When set alongside `--target-variants`, the k-mer index is queried directly for each undetected TI variant. Results appear with `call_source=RESCUED` (some evidence found below thresholds) or `call_source=NO_EVIDENCE` (zero supporting k-mers). Sub-threshold calls that match a TI target are marked `call_source=SUBTHRESHOLD`. This is a boolean flag: present means enabled.

**Example 1: Full tumour-informed monitoring with rescue**
```bash
kam run --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --target-variants patient_042_mutations.vcf \
  --ti-rescue \
  --output-format-override tsv,vcf
```

**Example 2: TI rescue with relaxed confidence**
```bash
kam run --r1 low_input_R1.fq.gz --r2 low_input_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --target-variants tumour_mutations.vcf \
  --ti-rescue \
  --min-confidence 0.95 \
  --min-alt-molecules 1
```

---

### SV and fusion options

#### `--sv-junctions`

FASTA of SV junction sequences to augment the k-mer allowlist. Required for detecting inversions and InvDel events whose breakpoint k-mers are absent from the panel targets. Not needed for large deletions, tandem duplications, or novel insertions (their k-mers come from the reference targets).

**Example 1: SV detection with junction sequences**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --sv-junctions sv_junctions.fa
```

**Example 2: SV detection alongside tumour-informed monitoring**
```bash
kam run --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --sv-junctions sv_junctions.fa \
  --target-variants patient_sv_mutations.vcf \
  --ti-position-tolerance-override 5
```

**Example 3: Combined SV and fusion detection**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --sv-junctions sv_junctions.fa \
  --fusion-targets fusion_targets.fa \
  --output-format-override tsv,vcf
```

---

#### `--fusion-targets`

FASTA of synthetic fusion target sequences for translocation/gene fusion detection. Each entry ID must follow the format `{name}__{chromA}:{startA}-{endA}__{chromB}:{startB}-{endB}__fusion`. K-mers from these sequences are added to the allowlist so fusion-spanning reads are captured. After normal target processing, fusion targets are walked and called separately using partner depth as the VAF denominator.

**Example 1: BCR-ABL fusion detection**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets twist_myeloid_panel.fa --output-dir results/ \
  --fusion-targets bcr_abl_fusions.fa
```

**Example 2: Multi-fusion panel**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --fusion-targets fusion_targets.fa \
  --sv-junctions sv_junctions.fa \
  --output-format-override tsv,vcf
```

**Example 3: Fusion detection with relaxed SV thresholds**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --fusion-targets fusion_targets.fa \
  --sv-min-confidence 0.90 \
  --sv-min-alt-molecules 1
```

---

#### `--junction-sequences`

FASTA of raw junction sequences for monitoring fusions or SVs observed in BAM/IGV. Any FASTA header format is accepted (no coordinate syntax required). Each sequence is added to the k-mer allowlist and walked as a standalone target with total library depth as the VAF denominator. Use this when you have the observed junction sequence but not exact genomic coordinates.

**Example 1: Monitor a known junction from IGV**
```bash
kam run --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --junction-sequences patient_042_observed_junctions.fa
```

**Example 2: Combine with SV junctions**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --sv-junctions sv_junctions.fa \
  --junction-sequences igv_junctions.fa
```

**Example 3: TI monitoring with junction sequence tracking**
```bash
kam run --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --target-variants tumour_mutations.vcf \
  --junction-sequences patient_fusion_junctions.fa \
  --ti-rescue
```

---

#### `--sv-min-confidence`

Minimum posterior probability for a structural variant PASS call. Applies to LargeDeletion, TandemDuplication, Inversion, InvDel, NovelInsertion, and Fusion types. Default: 0.95. Lower than the SNV default (0.99) because a large structural event with 2 supporting molecules is qualitatively stronger evidence than 2 molecules supporting a single-base change.

**Example 1: Default SV confidence**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --sv-junctions sv_junctions.fa \
  --sv-min-confidence 0.95
```

**Example 2: Relaxed for discovery**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --sv-junctions sv_junctions.fa \
  --sv-min-confidence 0.80
```

---

#### `--sv-min-alt-molecules`

Minimum alt-supporting molecules for a structural variant PASS call. Default: 1. A single-molecule SV call can be meaningful in monitoring mode where the target allele is pre-specified.

**Example 1: Single molecule sufficient (monitoring)**
```bash
kam run --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --sv-junctions sv_junctions.fa \
  --target-variants sv_mutations.vcf \
  --sv-min-alt-molecules 1
```

**Example 2: Require 2 molecules for discovery**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --sv-junctions sv_junctions.fa \
  --sv-min-alt-molecules 2
```

---

#### `--sv-strand-bias-threshold-override`

Fisher p-value threshold for strand bias on SV-type variants. Default: 1.0 (disabled). Inversion junction reads are structurally strand-biased due to directional path walking, so the standard threshold is inappropriate for SV paths. Set to a value below 1.0 to enable filtering, or 0.0 to apply the same threshold as SNVs/indels. Overrides the config file value.

**Example 1: Disabled (default)**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --sv-junctions sv_junctions.fa \
  --sv-strand-bias-threshold-override 1.0
```

**Example 2: Apply SNV-level strand bias to SVs**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --sv-junctions sv_junctions.fa \
  --sv-strand-bias-threshold-override 0.01
```

---

### ML options

#### `--ml-model`

Built-in model name for variant re-scoring. The model is bundled with the binary and works without any external files. Adds `ml_prob` and `ml_filter` columns to the output. Mutually exclusive with `--custom-ml-model`.

Available built-in models (use `kam models list` to see the current list):
- `single-strand-v1`
- `twist-duplex-v1`
- `twist-duplex-v2`

**Example 1: Duplex-trained model**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --ml-model twist-duplex-v2
```

**Example 2: Single-strand model for simplex data**
```bash
kam run --r1 simplex_R1.fq.gz --r2 simplex_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --ml-model single-strand-v1
```

**Example 3: ML re-scoring with TI monitoring**
```bash
kam run --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --target-variants tumour_mutations.vcf \
  --ml-model twist-duplex-v2 \
  --output-format-override tsv,vcf
```

---

#### `--custom-ml-model`

Path to a custom ONNX model file. The companion `.json` metadata file must exist at the same path (e.g. `my_model.onnx` requires `my_model.json`). Mutually exclusive with `--ml-model`.

**Example 1: Custom trained model**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --custom-ml-model /models/lab_trained_v3.onnx
```

**Example 2: Panel-specific model**
```bash
kam run --config configs/discovery.toml \
  --custom-ml-model /ref/models/pancancer_panel_rescorer.onnx
```

---

### Output options

#### `--output-format-override`

Comma-separated list of output formats. Supported values: `tsv`, `csv`, `json`, `vcf`. When multiple formats are given, each is written to `<output_base>.<ext>` in the output directory. Overrides the config file value.

**Example 1: TSV only (default)**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --output-format-override tsv
```

**Example 2: TSV and VCF together**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --output-format-override tsv,vcf
```

**Example 3: All four formats**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --output-format-override tsv,csv,json,vcf
```

---

#### `--qc-output`

Path for a merged QC JSON file. Per-stage QC JSONs (`assemble_qc.json`, `index_qc.json`, `pathfind_qc.json`, `call_qc.json`) are always written to the output directory regardless of this flag.

**Example 1: Specify merged QC path**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --qc-output results/merged_qc.json
```

**Example 2: QC output alongside results**
```bash
kam run --config configs/discovery.toml \
  --qc-output /data/qc/sample_042_qc.json
```

---

### Logging options

#### `--log-dir`

Directory for structured log output.

**Example 1: Local log directory**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --log-dir results/logs/
```

**Example 2: Central logging location**
```bash
kam run --config configs/discovery.toml \
  --log-dir /var/log/kam/run_2025_04_16/
```

---

#### `--log`

Enable specific log sinks. Can be repeated to enable multiple channels.

**Example 1: UMI-level logging**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --log-dir results/logs/ --log umi
```

**Example 2: Multiple log channels**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --log-dir results/logs/ --log umi --log family
```

---

### Performance options

#### `--threads`

Number of worker threads.

**Example 1: Use 8 threads**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --threads 8
```

**Example 2: Single thread for reproducibility debugging**
```bash
kam run --config configs/discovery.toml --threads 1
```

---

## kam assemble

Parse paired FASTQ files and produce a molecules bincode file. Writes `assemble_qc.json` to the same directory as `--output`.

```
kam assemble --r1 <R1.fastq.gz> --r2 <R2.fastq.gz> --output <molecules.bin> [OPTIONS]
```

### Required arguments

#### `--r1`

R1 input FASTQ file.

**Example 1: Gzip-compressed input**
```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin
```

**Example 2: Absolute path**
```bash
kam assemble --r1 /data/sequencing/run_001/sample_R1.fq.gz \
  --r2 /data/sequencing/run_001/sample_R2.fq.gz \
  --output /scratch/kam/sample_molecules.bin
```

---

#### `--r2`

R2 input FASTQ file. Must match `--r1` in read order and count.

**Example: Standard paired-end**
```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin
```

---

#### `--output`

Output file path for assembled molecules (bincode format).

**Example 1: Local output**
```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output results/molecules.bin
```

**Example 2: Scratch storage**
```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output /scratch/kam/batch_01/molecules.bin
```

---

### Optional arguments

#### `--chemistry`

Chemistry preset. Default: `twist-umi-duplex`.

**Example 1: Explicit Twist chemistry**
```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin --chemistry twist-umi-duplex
```

---

#### `--min-umi-quality`

Minimum Phred quality threshold for UMI bases. Default: 20. Read pairs with any UMI base below this are dropped.

**Example 1: Increase to Q30**
```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin --min-umi-quality 30
```

**Example 2: Disable quality filtering**
```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin --min-umi-quality 0
```

---

#### `--min-family-size`

Minimum number of reads per UMI family. Default: 1.

**Example 1: Require at least 2 reads**
```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin --min-family-size 2
```

**Example 2: Strict families for high-depth data**
```bash
kam assemble --r1 high_depth_R1.fq.gz --r2 high_depth_R2.fq.gz \
  --output molecules.bin --min-family-size 3
```

---

#### `--min-template-length`

Minimum template length in bases. Not set by default.

**Example: Filter short templates**
```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin --min-template-length 30
```

---

#### `--log-dir`

Directory for log output.

**Example: Enable logging**
```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin --log-dir logs/
```

---

#### `--log`

Enable specific log sinks. Repeatable.

**Example: Multiple log channels**
```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin --log-dir logs/ --log umi --log family
```

---

#### `--threads`

Number of worker threads.

**Example: Use 4 threads**
```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin --threads 4
```

---

## kam index

Build a k-mer index from assembled molecules against target sequences. Writes `index_qc.json` to the same directory as `--output`.

```
kam index --input <molecules.bin> --targets <targets.fa> --output <index.bin> [OPTIONS]
```

### Required arguments

#### `--input`

Consensus molecules file produced by `kam assemble`.

**Example: Standard indexing**
```bash
kam index --input molecules.bin --targets panel.fa --output index.bin
```

---

#### `--targets`

Target sequences FASTA file. The same file used in `kam run` or passed to `kam pathfind`.

**Example: Panel targets**
```bash
kam index --input molecules.bin --targets twist_cfdna_pancancer.fa \
  --output index.bin
```

---

#### `--output`

Output k-mer index file (bincode format).

**Example: Write to scratch**
```bash
kam index --input molecules.bin --targets panel.fa \
  --output /scratch/kam/index.bin
```

---

### Optional arguments

#### `-k` / `--kmer-size`

K-mer size. Must be 1 to 31. Default: 31.

**Example 1: Default k=31**
```bash
kam index --input molecules.bin --targets panel.fa \
  --output index.bin -k 31
```

**Example 2: Shorter k for small targets**
```bash
kam index --input molecules.bin --targets short_amplicons.fa \
  --output index.bin -k 21
```

**Example 3: Explicit long form**
```bash
kam index --input molecules.bin --targets panel.fa \
  --output index.bin --kmer-size 25
```

---

#### `--sv-junctions`

FASTA of SV junction sequences to augment the k-mer allowlist. Same purpose as the `kam run` flag.

**Example: Add SV junction k-mers**
```bash
kam index --input molecules.bin --targets panel.fa \
  --output index.bin --sv-junctions sv_junctions.fa
```

---

#### `--junction-sequences`

FASTA of raw junction sequences for augmenting the allowlist. Any header format is accepted.

**Example: Add observed junction sequences**
```bash
kam index --input molecules.bin --targets panel.fa \
  --output index.bin --junction-sequences patient_042_junctions.fa
```

---

## kam pathfind

Walk de Bruijn graph paths through the k-mer index for each target. Writes `pathfind_qc.json` to the same directory as `--output`.

```
kam pathfind --index <index.bin> --targets <targets.fa> --output <paths.bin> [OPTIONS]
```

### Required arguments

#### `--index`

K-mer index file produced by `kam index`.

**Example: Standard path walking**
```bash
kam pathfind --index index.bin --targets panel.fa --output paths.bin
```

---

#### `--targets`

Target sequences FASTA file. Must be the same file used during indexing.

**Example: Panel targets**
```bash
kam pathfind --index index.bin --targets twist_cfdna_pancancer.fa \
  --output paths.bin
```

---

#### `--output`

Output scored paths file (bincode format).

**Example: Write to results directory**
```bash
kam pathfind --index index.bin --targets panel.fa \
  --output results/paths.bin
```

---

### Optional arguments

#### `-k` / `--kmer-size`

K-mer size. Overrides the value inferred from the index. By default, k is inferred from the stored k-mer values. Use this flag when inference fails (e.g. all k-mers begin with A) or to confirm the expected k explicitly.

**Example 1: Explicit k to match indexing**
```bash
kam pathfind --index index.bin --targets panel.fa \
  --output paths.bin -k 31
```

**Example 2: Non-default k**
```bash
kam pathfind --index index.bin --targets panel.fa \
  --output paths.bin --kmer-size 25
```

---

## kam call

Call variants from scored paths. Writes `call_qc.json` to the same directory as `--output`.

```
kam call --paths <paths.bin> --output <variants.tsv> [OPTIONS]
```

### Required arguments

#### `--paths`

Scored paths file produced by `kam pathfind`.

**Example: Standard variant calling**
```bash
kam call --paths paths.bin --output variants.tsv
```

---

#### `--output`

Output variant calls file. The extension is informational; the actual format is controlled by `--output-format`.

**Example 1: TSV output**
```bash
kam call --paths paths.bin --output variants.tsv
```

**Example 2: VCF output**
```bash
kam call --paths paths.bin --output variants.vcf --output-format vcf
```

---

### Optional arguments

#### `--output-format`

Output format(s), comma-separated. Supported values: `tsv`, `csv`, `json`, `vcf`. When multiple formats are given, each is written to `<output_base>.<ext>` using the output file base name. Default: `tsv`.

**Example 1: TSV only**
```bash
kam call --paths paths.bin --output variants.tsv --output-format tsv
```

**Example 2: TSV and VCF**
```bash
kam call --paths paths.bin --output variants.tsv --output-format tsv,vcf
# Writes: variants.tsv and variants.vcf
```

**Example 3: All formats**
```bash
kam call --paths paths.bin --output variants.tsv --output-format tsv,csv,json,vcf
# Writes: variants.tsv, variants.csv, variants.json, variants.vcf
```

---

#### `--min-confidence`

Minimum posterior probability for a variant call. Default: 0.99 (from CallerConfig).

**Example 1: Default confidence**
```bash
kam call --paths paths.bin --output variants.tsv --min-confidence 0.99
```

**Example 2: Relaxed for discovery**
```bash
kam call --paths paths.bin --output variants.tsv --min-confidence 0.90
```

---

#### `--strand-bias-threshold`

Strand bias Fisher's exact p-value cutoff. Default: 0.01 (from CallerConfig).

**Example 1: Default threshold**
```bash
kam call --paths paths.bin --output variants.tsv --strand-bias-threshold 0.01
```

**Example 2: Stricter for FFPE**
```bash
kam call --paths paths.bin --output variants.tsv --strand-bias-threshold 0.05
```

---

#### `--min-alt-molecules`

Minimum alt-supporting molecules. Default: 2 (from CallerConfig).

**Example 1: Require 3 molecules**
```bash
kam call --paths paths.bin --output variants.tsv --min-alt-molecules 3
```

**Example 2: Allow single-molecule calls in monitoring**
```bash
kam call --paths paths.bin --output variants.tsv \
  --min-alt-molecules 1 --target-variants tumour_mutations.vcf
```

---

#### `--min-alt-duplex`

Minimum variant-specific duplex molecules. Default: 0 (disabled).

**Example 1: Require duplex confirmation**
```bash
kam call --paths paths.bin --output variants.tsv --min-alt-duplex 1
```

**Example 2: Disabled (default)**
```bash
kam call --paths paths.bin --output variants.tsv --min-alt-duplex 0
```

---

#### `--sv-min-confidence`

Minimum posterior probability for SV-type calls. Default: 0.95.

**Example: Relaxed SV confidence**
```bash
kam call --paths paths.bin --output variants.tsv --sv-min-confidence 0.85
```

---

#### `--sv-min-alt-molecules`

Minimum alt molecules for SV-type calls. Default: 1.

**Example: Require 2 molecules for SV calls**
```bash
kam call --paths paths.bin --output variants.tsv --sv-min-alt-molecules 2
```

---

#### `--sv-strand-bias-threshold`

Fisher p-value threshold for strand bias on SV-type variants. Default: 1.0 (disabled).

**Example 1: Disabled (default)**
```bash
kam call --paths paths.bin --output variants.tsv --sv-strand-bias-threshold 1.0
```

**Example 2: Enable SV strand bias filtering**
```bash
kam call --paths paths.bin --output variants.tsv --sv-strand-bias-threshold 0.01
```

---

#### `--max-vaf`

Maximum VAF for a PASS call. Calls above this are labelled `HighVaf`. Not set by default.

**Example: Exclude germline variants**
```bash
kam call --paths paths.bin --output variants.tsv --max-vaf 0.35
```

---

#### `--target-variants`

VCF file of expected somatic variants for tumour-informed monitoring. Same behaviour as the `kam run` flag.

**Example: Tumour-informed calling**
```bash
kam call --paths paths.bin --output variants.tsv \
  --target-variants tumour_biopsy.vcf --max-vaf 0.35
```

---

#### `--ti-position-tolerance`

Position tolerance in bp for tumour-informed matching. Default: 0.

**Example: Allow 5 bp tolerance for SVs**
```bash
kam call --paths paths.bin --output variants.tsv \
  --target-variants sv_mutations.vcf --ti-position-tolerance 5
```

---

#### `--ml-model`

Built-in model name for variant re-scoring. Same behaviour as the `kam run` flag.

**Example: Apply ML re-scoring**
```bash
kam call --paths paths.bin --output variants.tsv \
  --ml-model twist-duplex-v2 --output-format tsv,vcf
```

---

#### `--custom-ml-model`

Path to a custom ONNX model file. Same behaviour as the `kam run` flag.

**Example: Custom model re-scoring**
```bash
kam call --paths paths.bin --output variants.tsv \
  --custom-ml-model /models/lab_trained_v3.onnx
```

---

## kam explore

Open an interactive REPL for inspecting bincode data files produced by kam. The file type (molecules, k-mer index, or variant calls) is detected automatically from the bincode header.

```
kam explore <FILE>
```

### Arguments

#### `<FILE>`

Path to a bincode file. Accepts any `.bin` file produced by `kam assemble`, `kam index`, `kam call`, or `kam run`.

### Commands

Once inside the explorer, type `help` to see all available commands. Key commands: `summary`, `head N`, `show <id>`, `filter <expression>`, `stats <field>`, `histogram <field>`, `export`, `quit`.

See the [Interactive explorer guide](guides/interactive-explorer.md) for full command reference and filter expression syntax.

**Example 1: Explore assembled molecules**
```bash
kam explore results/molecules.bin
```

```
kam explore v0.3.0
Loaded: molecules.bin (612,847 molecules, 4,821,044 input reads)
Type 'help' for commands, 'quit' to exit.

molecules> summary
  Total molecules: 612,847
  Duplex families: 489,478 (79.8%)
  ...
```

**Example 2: Explore variant calls and filter**
```bash
kam explore results/variants.bin
```

```
variants> filter filter == PASS AND vaf > 0.001
  Matched 12 variants

variants> export --format tsv --output pass_calls.tsv
  Exported 12 variants to pass_calls.tsv
```

**Example 3: Explore a k-mer index**
```bash
kam explore results/index.bin
```

```
kmers> summary
  Total k-mers:   1,247,891
  K-mer size:     31
  Target k-mers:  89,412
  Junction k-mers: 1,247
  ...
```

---

## kam models list

List all built-in ML models available for variant re-scoring.

```
kam models list
```

**Example output:**
```
Built-in models:
  single-strand-v1
  twist-duplex-v1
  twist-duplex-v2
```

Use these names with `--ml-model` on `kam run` or `kam call`.

---

## Output format detail

### Multiple formats

Passing a comma-separated list to `--output-format` (or `--output-format-override` for `run`) writes each format to a separate file. The output path base (without extension) is used:

```bash
kam call --paths paths.bin --output variants.tsv --output-format tsv,vcf
# Writes: variants.tsv and variants.vcf
```

### VCF coordinate placement

When a target ID has the format `chrN:START-END` (0-based, half-open interval), the VCF writer places variants at their correct genomic coordinates:

- CHROM = `chrN`
- POS = `START + offset` (1-based; the offset within the full target sequence where the allele differs)
- REF and ALT = minimal left-normalised alleles

When the target ID cannot be parsed as a coordinate, the writer falls back to:
- CHROM = target ID as-is
- POS = 1

---

## Environment variables

None currently used.

---

## Full pipeline examples

### Discovery mode

Standard variant discovery across a gene panel with no prior knowledge of expected variants.

```bash
kam run \
  --r1 patient_042_R1.fq.gz \
  --r2 patient_042_R2.fq.gz \
  --targets twist_cfdna_pancancer.fa \
  --output-dir results/patient_042/ \
  --output-format-override tsv,vcf
```

### Tumour-informed monitoring mode

Track known somatic mutations from a prior tissue biopsy across serial plasma time points.

```bash
kam run \
  --r1 plasma_tp3_R1.fq.gz \
  --r2 plasma_tp3_R2.fq.gz \
  --targets twist_cfdna_pancancer.fa \
  --output-dir results/patient_042_tp3/ \
  --target-variants patient_042_mutations.vcf \
  --ti-rescue \
  --max-vaf 0.35 \
  --output-format-override tsv,vcf
```

### TI monitoring with SV tracking

Monitor SNVs, indels, and structural variants from the original tumour biopsy.

```bash
kam run \
  --r1 plasma_R1.fq.gz \
  --r2 plasma_R2.fq.gz \
  --targets panel.fa \
  --output-dir results/ \
  --target-variants patient_sv_mutations.vcf \
  --sv-junctions sv_junctions.fa \
  --junction-sequences patient_observed_junctions.fa \
  --ti-rescue \
  --ti-position-tolerance-override 5 \
  --max-vaf 0.35 \
  --output-format-override tsv,vcf
```

### Run with duplex confirmation required

High-specificity mode requiring at least one duplex molecule at the variant site.

```bash
kam run \
  --r1 high_depth_R1.fq.gz \
  --r2 high_depth_R2.fq.gz \
  --targets panel.fa \
  --output-dir results/ \
  --min-alt-duplex 1 \
  --output-format-override tsv,vcf
```

### Config file with CLI overrides

Use a config file for baseline settings and override specific values per sample.

```bash
kam run \
  --config configs/tumour-informed.toml \
  --r1 new_plasma_R1.fq.gz \
  --r2 new_plasma_R2.fq.gz \
  --target-variants new_patient_mutations.vcf \
  --output-dir results/new_patient/
```

### Fusion detection

Detect gene fusions alongside standard variants.

```bash
kam run \
  --r1 sample_R1.fq.gz \
  --r2 sample_R2.fq.gz \
  --targets twist_myeloid_panel.fa \
  --output-dir results/ \
  --fusion-targets bcr_abl_fusions.fa \
  --sv-min-confidence 0.90 \
  --output-format-override tsv,vcf
```

### Stage-by-stage pipeline

Run each stage independently for debugging or Nextflow integration.

```bash
# 1. Assemble molecules
kam assemble \
  --r1 sample_R1.fq.gz \
  --r2 sample_R2.fq.gz \
  --output molecules.bin \
  --min-family-size 2

# 2. Build k-mer index
kam index \
  --input molecules.bin \
  --targets panel.fa \
  --output index.bin \
  -k 31 \
  --sv-junctions sv_junctions.fa

# 3. Walk graph paths
kam pathfind \
  --index index.bin \
  --targets panel.fa \
  --output paths.bin

# 4. Call variants
kam call \
  --paths paths.bin \
  --output variants.tsv \
  --output-format tsv,vcf \
  --min-confidence 0.99 \
  --max-vaf 0.35 \
  --target-variants tumour_biopsy.vcf
```

### ML re-scoring

Apply a built-in gradient boosting model for variant re-scoring.

```bash
kam run \
  --r1 sample_R1.fq.gz \
  --r2 sample_R2.fq.gz \
  --targets panel.fa \
  --output-dir results/ \
  --ml-model twist-duplex-v2 \
  --output-format-override tsv,vcf
```

---

## Notes

### File formats

All intermediate binary files (`*.bin`) are `bincode`-serialised with a header identifying the file type and the kam version. Passing a file produced by one version of kam to a different version will fail with an error if the serialisation format changed.

### Determinism

All pipeline stages produce deterministic output for a given input, regardless of the number of threads. Internal HashMap ordering is not relied upon for output: molecules are sorted by UMI hash, variants are sorted by target ID before writing.

### Error messages

All errors print to stderr. The exit code is 1 on any error, 0 on success. Stage progress is reported to stderr with `[stage]` prefixes, e.g. `[assemble] molecules=351341 duplex=66552`.

### Flag naming conventions

`kam run` uses `_override` suffixes on certain flags to distinguish CLI overrides from config file values. The underlying behaviour is identical:

| `kam run` flag | Standalone subcommand equivalent |
|---|---|
| `--chemistry-override` | `--chemistry` (assemble) |
| `--min-umi-quality-override` | `--min-umi-quality` (assemble) |
| `--min-family-size-override` | `--min-family-size` (assemble) |
| `--kmer-size-override` / `-k` | `--kmer-size` / `-k` (index) |
| `--output-format-override` | `--output-format` (call) |
| `--ti-position-tolerance-override` | `--ti-position-tolerance` (call) |
| `--sv-strand-bias-threshold-override` | `--sv-strand-bias-threshold` (call) |
