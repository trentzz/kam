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

## Input file formats

Every file-path flag in kam expects a specific format. This section shows what each input file looks like so you can prepare your own.

### Targets FASTA (`--targets`)

Each entry defines one target window. The header encodes a genomic coordinate in `chrN:START-END` format (0-based, half-open). This enables correct VCF coordinate placement. Sequences are typically 100 bp of reference extracted from the genome.

```
>chr17:7674220-7674320
ATGCAGTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr7:55241607-55241707
TGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAAT
GCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGC
>BRCA1_exon10:chr17:43091434-43091534
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
```

Headers can optionally include a descriptive prefix before the coordinate (e.g. `BRCA1_exon10:chr17:...`). The coordinate portion is parsed for VCF output. If no coordinate is parseable, the full header is used as CHROM and POS is set to 1.

### SV junctions FASTA (`--sv-junctions`)

Each entry is a synthetic sequence spanning a structural variant breakpoint. The header encodes the coordinate, SV type, and event size. The sequence concatenates left-flank and right-flank reference bases around the breakpoint.

```
>chr17:7674100-7674300_DEL_100bp
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
>chr7:55241500-55241800_INV_200bp
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
>chr12:25380100-25380400_DUP_150bp
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
```

Required for detecting inversions and InvDel events whose breakpoint k-mers are absent from the panel targets. Not needed for large deletions, tandem duplications, or novel insertions (their k-mers come from the reference targets). Generate this file using varforge or a BED of known breakpoints.

### Fusion targets FASTA (`--fusion-targets`)

Each entry is a synthetic sequence spanning a gene fusion breakpoint. The header **must** follow this exact format:

```
>{name}__{chromA}:{startA}-{endA}__{chromB}:{startB}-{endB}__fusion
```

The sequence concatenates partner A's breakpoint-adjacent segment and partner B's breakpoint-adjacent segment.

```
>BCR_ABL1__chr22:23632500-23632550__chr9:130854000-130854050__fusion
ACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
TGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATG
>EML4_ALK__chr2:42396490-42396540__chr2:29416089-29416139__fusion
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
```

Fields: `name` is a human-readable fusion name. `chromA:startA-endA` and `chromB:startB-endB` are the genomic coordinates of each partner segment. The `__fusion` suffix marks the entry as a fusion target.

### Junction sequences FASTA (`--junction-sequences`)

Raw junction sequences copied from BAM or IGV. Any header format is accepted. Each sequence is added to the k-mer allowlist and walked as a standalone target. Use this when you have the observed junction sequence but not exact genomic coordinates.

```
>fusion_observed_in_igv
ACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
TGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATG
>sv_breakpoint_from_bam
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
```

### Target variants VCF (`--target-variants`)

A standard VCF (v4.3) listing expected somatic variants from a prior tissue biopsy. Only the `#CHROM`, `POS`, `REF`, and `ALT` columns are used for matching. `QUAL`, `FILTER`, and `INFO` can be `.` (empty).

```
##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr17	7674220	.	C	T	.	.	.
chr7	55241670	.	ACGT	A	.	.	.
chr12	25380275	.	G	A	.	.	.
```

### Config file (`--config`)

A TOML file. CLI flags always override config file values, which override built-in defaults. See the [Configuration Reference](configuration.md) for the full schema.

```toml
[input]
r1 = "sample_R1.fq.gz"
r2 = "sample_R2.fq.gz"
targets = "panel.fa"

[output]
output_dir = "results/"
output_format = "tsv,vcf"

[chemistry]
umi_length = 5
skip_length = 2
duplex = true
min_umi_quality = 20

[calling]
min_confidence = 0.99
strand_bias_threshold = 0.01
min_alt_molecules = 2
max_vaf = 0.35
```

### Custom ML model (`--custom-ml-model`)

An ONNX model file. A companion `.json` metadata file must exist at the same path (e.g. `my_model.onnx` requires `my_model.json`).

The metadata JSON contains:

```json
{
  "version": "lab-trained-v3",
  "feature_names": [
    "vaf", "nref", "nalt", "ndupalt", "nsimalt", "sbp", "conf",
    "ref_len", "alt_len", "duplex_frac", "has_duplex", "ci_width",
    "alt_depth", "log_nalt", "log_nref", "log_alt_depth", "log_vaf",
    "vaf_times_conf", "vaf_times_nalt", "nalt_over_conf",
    "ci_width_rel", "snr", "conf_sq", "nalt_sq", "vaf_sq",
    "ref_alt_len_ratio", "indel_size", "duplex_enrichment",
    "simplex_only_frac", "conf_above_99", "conf_above_999",
    "sbp_above_05", "variant_class_enc"
  ],
  "ml_pass_threshold": 0.5,
  "variant_class_map": {
    "SNV": 0,
    "Insertion": 1,
    "Deletion": 2,
    "MNV": 3,
    "Complex": 4,
    "LargeDeletion": 5,
    "TandemDuplication": 6,
    "Inversion": 7,
    "Fusion": 8,
    "InvDel": 9,
    "NovelInsertion": 10
  }
}
```

`feature_names` lists the features in the exact order the model expects. `ml_pass_threshold` is the probability above which a call is labelled `ML_PASS`. `variant_class_map` maps variant type strings to integer encodings.

---

## kam run

Run the complete pipeline: assemble, index, pathfind, call.

```
kam run --r1 <R1.fq.gz> --r2 <R2.fq.gz> --targets <targets.fa> --output-dir <DIR> [OPTIONS]
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

```bash
kam run --r1 patient_042_R1.fq.gz --r2 patient_042_R2.fq.gz \
  --targets panel.fa --output-dir results/
```

---

#### `--r2`

R2 input FASTQ file. Must be paired with `--r1`. Same format requirements.

```bash
kam run --r1 ctdna_R1.fq.gz --r2 ctdna_R2.fq.gz \
  --targets panel.fa --output-dir results/
```

---

#### `--targets`

Target sequences FASTA file. Each entry defines one target window. See [Targets FASTA](#targets-fasta---targets) above for the file format.

```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets twist_cfdna_pancancer.fa --output-dir results/
```

---

#### `--output-dir`

Directory for all pipeline outputs. Created if it does not exist. Stage-specific QC JSON files, intermediate binaries, and final variant calls are all written here.

```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir /scratch/kam_results/patient_042/
```

---

#### `--config`

Path to a TOML configuration file. When provided, pipeline parameters are loaded from the file first. CLI flags then override individual values. See [Config file](#config-file---config) above for the format, and the [Configuration Reference](configuration.md) for the full schema.

**Config file only:**
```bash
kam run --config configs/discovery.toml
```

**Config with CLI overrides:**
```bash
kam run --config configs/tumour-informed.toml \
  --r1 new_sample_R1.fq.gz --r2 new_sample_R2.fq.gz \
  --target-variants patient_042_mutations.vcf \
  --output-dir results/new_patient/
```

---

### Assembly options

#### `--chemistry-override`

Chemistry preset name. Overrides the config file value. Currently only `twist-umi-duplex` is supported, which sets UMI length = 5 bp, skip length = 2 bp, and template start = 7 bp.

```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --chemistry-override twist-umi-duplex
```

---

#### `--min-umi-quality-override`

Minimum Phred quality for UMI bases. Read pairs with any UMI base below this quality are dropped before grouping. Default: 20. Set to 0 to disable quality filtering. Overrides the config file value.

**Increase stringency to Q30:**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-umi-quality-override 30
```

**Disable UMI quality filtering:**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-umi-quality-override 0
```

---

#### `--min-family-size-override`

Minimum reads per UMI family. Families with fewer reads are discarded. Default: 1. Setting to 2 eliminates true singletons but reduces sensitivity at low depth. Overrides the config file value.

```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-family-size-override 2
```

---

#### `--min-template-length`

Minimum template length in bases. Templates shorter than this on either R1 or R2 are dropped. Not set by default. Useful for filtering adapter dimers and very short inserts.

```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-template-length 30
```

---

### Indexing options

#### `-k` / `--kmer-size-override`

K-mer length for indexing. Must be 1 to 31. Larger k improves specificity (fewer spurious paths) but requires reads that fully span the k-mer window. At k=31, reads shorter than 31 bp produce no k-mers. Default: 31. Overrides the config file value.

**Default k=31 for standard panels:**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ -k 31
```

**Shorter k for short target windows or low-depth samples:**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets short_amplicons.fa --output-dir results/ -k 21
```

---

### Calling options

#### `--min-confidence`

Minimum posterior probability for a PASS call. Derived from a binomial likelihood ratio between the observed VAF and the background error rate. Calls below this threshold are labelled `LowConfidence`. Default: 0.99.

**Relaxed for discovery in low-depth samples:**
```bash
kam run --r1 low_depth_R1.fq.gz --r2 low_depth_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-confidence 0.90
```

**Ultra-stringent for clinical reporting:**
```bash
kam run --config configs/tumour-informed.toml \
  --min-confidence 0.999
```

---

#### `--strand-bias-threshold`

Fisher's exact test p-value cutoff for strand bias. Calls with p below this threshold are labelled `StrandBias`. Genuine somatic variants appear on both strands. Artefacts (oxidative damage, end-repair) tend to be strand-specific. Default: 0.01.

**Stricter filter for FFPE samples (more artefacts):**
```bash
kam run --r1 ffpe_R1.fq.gz --r2 ffpe_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --strand-bias-threshold 0.05
```

---

#### `--min-alt-molecules`

Minimum molecules supporting the alt allele for a PASS call. A single-molecule call (n_alt=1) is accepted only if duplex-confirmed. Default: 2.

**High specificity for clinical reporting:**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-alt-molecules 3
```

**Allow single-molecule calls in monitoring mode:**
```bash
kam run --r1 ctdna_R1.fq.gz --r2 ctdna_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --target-variants tumour_mutations.vcf \
  --min-alt-molecules 1
```

---

#### `--min-alt-duplex`

Minimum variant-specific duplex molecules. Counts only duplex molecules whose k-mers overlap the variant site (alt-path-specific k-mers), not flanking anchor k-mers. Default: 0 (disabled). Set to 1 to require at least one duplex confirmation at the variant site. Recommended only at duplex fractions above 15% or depth above 5M reads.

```bash
kam run --r1 high_depth_R1.fq.gz --r2 high_depth_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --min-alt-duplex 1
```

---

#### `--max-vaf`

Maximum VAF for a PASS call. Calls with VAF above this are labelled `HighVaf`. Not set by default. For ctDNA somatic calling, germline heterozygous variants have VAF around 0.5. Setting `--max-vaf 0.35` eliminates them while retaining all low-VAF somatic calls.

```bash
kam run --r1 ctdna_R1.fq.gz --r2 ctdna_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --max-vaf 0.35
```

---

#### `--target-variants`

VCF file of expected somatic variants for tumour-informed monitoring mode. When set, only calls whose `(CHROM, POS, REF, ALT)` matches an entry in this VCF are marked PASS. All other quality-passing calls are labelled `NotTargeted`. See [Target variants VCF](#target-variants-vcf---target-variants) above for the file format.

**Standard tumour-informed monitoring:**
```bash
kam run --r1 plasma_tp2_R1.fq.gz --r2 plasma_tp2_R2.fq.gz \
  --targets twist_cfdna_pancancer.fa --output-dir results/tp2/ \
  --target-variants tumour_biopsy_snvs.vcf \
  --max-vaf 0.35
```

**Combined with TI rescue for complete tracking:**
```bash
kam run --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --target-variants patient_042_mutations.vcf \
  --ti-rescue \
  --output-format-override tsv,vcf
```

---

#### `--ti-position-tolerance-override`

Position tolerance in base pairs for tumour-informed matching. Default: 0 (exact coordinate matching). When greater than 0, a call also passes the TI filter if its position is within this many bp of any target variant position, regardless of REF/ALT. Useful for large SVs with breakpoint ambiguity.

**Allow 1 bp tolerance for indel left-alignment differences:**
```bash
kam run --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --target-variants tumour_mutations.vcf \
  --ti-position-tolerance-override 1
```

**Large tolerance for SV breakpoints:**
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

```bash
kam run --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --target-variants patient_042_mutations.vcf \
  --ti-rescue \
  --min-confidence 0.95 \
  --min-alt-molecules 1 \
  --output-format-override tsv,vcf
```

---

### SV and fusion options

#### `--sv-junctions`

FASTA of SV junction sequences to augment the k-mer allowlist. See [SV junctions FASTA](#sv-junctions-fasta---sv-junctions) above for the file format.

Required for detecting inversions and InvDel events. Not needed for large deletions, tandem duplications, or novel insertions (their k-mers come from the reference targets).

**SV detection with junction sequences:**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --sv-junctions sv_junctions.fa
```

**Combined SV and fusion detection:**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --sv-junctions sv_junctions.fa \
  --fusion-targets fusion_targets.fa \
  --output-format-override tsv,vcf
```

---

#### `--fusion-targets`

FASTA of synthetic fusion target sequences for translocation/gene fusion detection. See [Fusion targets FASTA](#fusion-targets-fasta---fusion-targets) above for the required header format.

K-mers from these sequences are added to the allowlist so fusion-spanning reads are captured. Fusion targets are walked and called separately using partner depth as the VAF denominator.

```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets twist_myeloid_panel.fa --output-dir results/ \
  --fusion-targets bcr_abl_fusions.fa \
  --sv-min-confidence 0.90
```

---

#### `--junction-sequences`

FASTA of raw junction sequences for monitoring fusions or SVs observed in BAM/IGV. See [Junction sequences FASTA](#junction-sequences-fasta---junction-sequences) above for the file format. Any FASTA header format is accepted.

Each sequence is added to the k-mer allowlist and walked as a standalone target with total library depth as the VAF denominator.

```bash
kam run --r1 plasma_R1.fq.gz --r2 plasma_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --junction-sequences patient_042_observed_junctions.fa \
  --target-variants tumour_mutations.vcf \
  --ti-rescue
```

---

#### `--sv-min-confidence`

Minimum posterior probability for a structural variant PASS call. Applies to LargeDeletion, TandemDuplication, Inversion, InvDel, NovelInsertion, and Fusion types. Default: 0.95. Lower than the SNV default (0.99) because a large structural event with 2 supporting molecules is qualitatively stronger evidence than 2 molecules supporting a single-base change.

```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --sv-junctions sv_junctions.fa \
  --sv-min-confidence 0.80
```

---

#### `--sv-min-alt-molecules`

Minimum alt-supporting molecules for a structural variant PASS call. Default: 1. A single-molecule SV call can be meaningful in monitoring mode where the target allele is pre-specified.

```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --sv-junctions sv_junctions.fa \
  --sv-min-alt-molecules 2
```

---

#### `--sv-strand-bias-threshold-override`

Fisher p-value threshold for strand bias on SV-type variants. Default: 1.0 (disabled). Inversion junction reads are structurally strand-biased due to directional path walking, so the standard threshold is inappropriate for SV paths. Set to a value below 1.0 to enable filtering. Overrides the config file value.

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

**Duplex-trained model:**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --ml-model twist-duplex-v2
```

**Single-strand model for simplex data:**
```bash
kam run --r1 simplex_R1.fq.gz --r2 simplex_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --ml-model single-strand-v1
```

---

#### `--custom-ml-model`

Path to a custom ONNX model file. The companion `.json` metadata file must exist at the same path (e.g. `my_model.onnx` requires `my_model.json`). See [Custom ML model](#custom-ml-model---custom-ml-model) above for the metadata format. Mutually exclusive with `--ml-model`.

```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --custom-ml-model /models/lab_trained_v3.onnx
```

---

### Output options

#### `--output-format-override`

Comma-separated list of output formats. Supported values: `tsv`, `csv`, `json`, `vcf`. When multiple formats are given, each is written to `<output_base>.<ext>` in the output directory. Overrides the config file value.

**TSV and VCF together:**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --output-format-override tsv,vcf
```

**All four formats:**
```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --output-format-override tsv,csv,json,vcf
```

---

#### `--qc-output`

Path for a merged QC JSON file. Per-stage QC JSONs (`assemble_qc.json`, `index_qc.json`, `pathfind_qc.json`, `call_qc.json`) are always written to the output directory regardless of this flag.

```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --qc-output results/merged_qc.json
```

---

### Logging options

#### `--log-dir`

Directory for structured log output.

```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --log-dir results/logs/
```

---

#### `--log`

Enable specific log sinks. Can be repeated to enable multiple channels.

```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --log-dir results/logs/ --log umi --log family
```

---

### Performance options

#### `--threads`

Number of worker threads.

```bash
kam run --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --targets panel.fa --output-dir results/ \
  --threads 8
```

---

## kam assemble

Parse paired FASTQ files and produce a molecules bincode file. Writes `assemble_qc.json` to the same directory as `--output`.

```
kam assemble --r1 <R1.fq.gz> --r2 <R2.fq.gz> --output <molecules.bin> [OPTIONS]
```

### Required arguments

#### `--r1`

R1 input FASTQ file.

#### `--r2`

R2 input FASTQ file. Must match `--r1` in read order and count.

#### `--output`

Output file path for assembled molecules (bincode format).

```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output results/molecules.bin
```

---

### Optional arguments

#### `--chemistry`

Chemistry preset. Default: `twist-umi-duplex`.

```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin --chemistry twist-umi-duplex
```

---

#### `--min-umi-quality`

Minimum Phred quality threshold for UMI bases. Default: 20. Read pairs with any UMI base below this are dropped.

```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin --min-umi-quality 30
```

---

#### `--min-family-size`

Minimum number of reads per UMI family. Default: 1.

```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin --min-family-size 2
```

---

#### `--min-template-length`

Minimum template length in bases. Not set by default.

```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin --min-template-length 30
```

---

#### `--log-dir`

Directory for log output.

```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin --log-dir logs/
```

---

#### `--log`

Enable specific log sinks. Repeatable.

```bash
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin --log-dir logs/ --log umi --log family
```

---

#### `--threads`

Number of worker threads.

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

#### `--targets`

Target sequences FASTA file. The same file used in `kam run` or passed to `kam pathfind`. See [Targets FASTA](#targets-fasta---targets) for the format.

#### `--output`

Output k-mer index file (bincode format).

```bash
kam index --input molecules.bin --targets panel.fa --output index.bin
```

---

### Optional arguments

#### `-k` / `--kmer-size`

K-mer size. Must be 1 to 31. Default: 31.

```bash
kam index --input molecules.bin --targets panel.fa \
  --output index.bin -k 25
```

---

#### `--sv-junctions`

FASTA of SV junction sequences to augment the k-mer allowlist. Same purpose as the `kam run` flag. See [SV junctions FASTA](#sv-junctions-fasta---sv-junctions) for the format.

```bash
kam index --input molecules.bin --targets panel.fa \
  --output index.bin --sv-junctions sv_junctions.fa
```

---

#### `--junction-sequences`

FASTA of raw junction sequences for augmenting the allowlist. Any header format is accepted. See [Junction sequences FASTA](#junction-sequences-fasta---junction-sequences) for the format.

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

#### `--targets`

Target sequences FASTA file. Must be the same file used during indexing.

#### `--output`

Output scored paths file (bincode format).

```bash
kam pathfind --index index.bin --targets panel.fa --output paths.bin
```

---

### Optional arguments

#### `-k` / `--kmer-size`

K-mer size. Overrides the value inferred from the index. By default, k is inferred from the stored k-mer values. Use this flag when inference fails (e.g. all k-mers begin with A) or to confirm the expected k explicitly.

```bash
kam pathfind --index index.bin --targets panel.fa \
  --output paths.bin -k 31
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

#### `--output`

Output variant calls file. The extension is informational; the actual format is controlled by `--output-format`.

```bash
kam call --paths paths.bin --output variants.tsv
```

---

### Optional arguments

#### `--output-format`

Output format(s), comma-separated. Supported values: `tsv`, `csv`, `json`, `vcf`. When multiple formats are given, each is written to `<output_base>.<ext>` using the output file base name. Default: `tsv`.

```bash
kam call --paths paths.bin --output variants.tsv --output-format tsv,vcf
# Writes: variants.tsv and variants.vcf
```

---

#### `--min-confidence`

Minimum posterior probability for a variant call. Default: 0.99 (from CallerConfig).

```bash
kam call --paths paths.bin --output variants.tsv --min-confidence 0.90
```

---

#### `--strand-bias-threshold`

Strand bias Fisher's exact p-value cutoff. Default: 0.01 (from CallerConfig).

```bash
kam call --paths paths.bin --output variants.tsv --strand-bias-threshold 0.05
```

---

#### `--min-alt-molecules`

Minimum alt-supporting molecules. Default: 2 (from CallerConfig).

```bash
kam call --paths paths.bin --output variants.tsv --min-alt-molecules 3
```

---

#### `--min-alt-duplex`

Minimum variant-specific duplex molecules. Default: 0 (disabled).

```bash
kam call --paths paths.bin --output variants.tsv --min-alt-duplex 1
```

---

#### `--sv-min-confidence`

Minimum posterior probability for SV-type calls. Default: 0.95.

```bash
kam call --paths paths.bin --output variants.tsv --sv-min-confidence 0.85
```

---

#### `--sv-min-alt-molecules`

Minimum alt molecules for SV-type calls. Default: 1.

```bash
kam call --paths paths.bin --output variants.tsv --sv-min-alt-molecules 2
```

---

#### `--sv-strand-bias-threshold`

Fisher p-value threshold for strand bias on SV-type variants. Default: 1.0 (disabled).

```bash
kam call --paths paths.bin --output variants.tsv --sv-strand-bias-threshold 0.01
```

---

#### `--max-vaf`

Maximum VAF for a PASS call. Calls above this are labelled `HighVaf`. Not set by default.

```bash
kam call --paths paths.bin --output variants.tsv --max-vaf 0.35
```

---

#### `--target-variants`

VCF file of expected somatic variants for tumour-informed monitoring. Same behaviour as the `kam run` flag. See [Target variants VCF](#target-variants-vcf---target-variants) for the format.

```bash
kam call --paths paths.bin --output variants.tsv \
  --target-variants tumour_biopsy.vcf --max-vaf 0.35
```

---

#### `--ti-position-tolerance`

Position tolerance in bp for tumour-informed matching. Default: 0.

```bash
kam call --paths paths.bin --output variants.tsv \
  --target-variants sv_mutations.vcf --ti-position-tolerance 5
```

---

#### `--ml-model`

Built-in model name for variant re-scoring. Same behaviour as the `kam run` flag.

```bash
kam call --paths paths.bin --output variants.tsv \
  --ml-model twist-duplex-v2 --output-format tsv,vcf
```

---

#### `--custom-ml-model`

Path to a custom ONNX model file. Same behaviour as the `kam run` flag. See [Custom ML model](#custom-ml-model---custom-ml-model) for the metadata format.

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

Each example below shows what input files you need, their content format, the exact command, and what output to expect.

### Discovery mode

Standard variant discovery across a gene panel with no prior knowledge of expected variants.

**Input files:**

1. Paired FASTQ files from sequencing
2. Targets FASTA with one entry per target window:

```
>chr17:7674220-7674320
ATGCAGTCGATCGATCGATCG...100bp reference...
>chr7:55241607-55241707
TGCAATGCAATGCAATGCAAT...100bp reference...
>TP53_exon7:chr17:7674100-7674200
GCTAGCTAGCTAGCTAGCTAG...100bp reference...
```

**Command:**
```bash
kam run \
  --r1 patient_042_R1.fq.gz \
  --r2 patient_042_R2.fq.gz \
  --targets twist_cfdna_pancancer.fa \
  --output-dir results/patient_042/ \
  --output-format-override tsv,vcf
```

**Output:** `results/patient_042/variants.tsv`, `results/patient_042/variants.vcf`, plus QC JSONs for each stage (`assemble_qc.json`, `index_qc.json`, `pathfind_qc.json`, `call_qc.json`).

---

### Tumour-informed monitoring mode

Track known somatic mutations from a prior tissue biopsy across serial plasma time points.

**Input files:**

1. Paired FASTQ files from plasma ctDNA sequencing
2. Targets FASTA (same panel as discovery)
3. Target variants VCF from the tumour biopsy:

```
##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr17	7674220	.	C	T	.	.	.
chr7	55241670	.	ACGT	A	.	.	.
chr12	25380275	.	G	A	.	.	.
```

**Command:**
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

**Output:** Only calls matching the truth VCF are marked PASS. All other calls are labelled `NotTargeted`. Rescued variants appear with `call_source=RESCUED` or `call_source=NO_EVIDENCE`.

---

### SV and fusion detection

Detect structural variants and gene fusions alongside standard SNVs and indels.

**Input files:**

1. Paired FASTQ files
2. Targets FASTA (standard panel targets)
3. SV junctions FASTA:

```
>chr17:7674100-7674300_DEL_100bp
ACGTACGTACGTACGTACGTACGT...50bp left flank...
TGCATGCATGCATGCATGCATGCA...50bp right flank...
>chr7:55241500-55241800_INV_200bp
GATCGATCGATCGATCGATCGATC...50bp left flank...
CTAGCTAGCTAGCTAGCTAGCTAGC...50bp right flank (reverse complement)...
```

4. Fusion targets FASTA:

```
>BCR_ABL1__chr22:23632500-23632550__chr9:130854000-130854050__fusion
ACGATCGATCG...50bp from BCR...TGCAATGCCAA...50bp from ABL1...
```

**Command:**
```bash
kam run \
  --r1 sample_R1.fq.gz \
  --r2 sample_R2.fq.gz \
  --targets twist_myeloid_panel.fa \
  --output-dir results/ \
  --sv-junctions sv_junctions.fa \
  --fusion-targets bcr_abl_fusions.fa \
  --sv-min-confidence 0.90 \
  --output-format-override tsv,vcf
```

**Output:** SNV/indel calls alongside SV and fusion calls in the same output files. SV calls use relaxed thresholds.

---

### TI monitoring with SV tracking

Monitor SNVs, indels, and structural variants from the original tumour biopsy, including observed junction sequences.

**Input files:**

1. Paired FASTQ files from plasma
2. Targets FASTA
3. Target variants VCF (SNVs, indels, and SV coordinates)
4. SV junctions FASTA
5. Junction sequences FASTA (from BAM/IGV):

```
>patient_042_fusion_observed
ACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
TGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATG
```

**Command:**
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

**Output:** Calls for SNVs/indels, SVs, and junction sequences. Only variants matching the truth VCF (within 5 bp tolerance) are marked PASS. Rescued variants report remaining evidence.

---

### Duplex confirmation required

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

---

### Config file with CLI overrides

Use a config file for baseline settings and override specific values per sample. See [Config file](#config-file---config) above for the format.

**Config file (`configs/tumour-informed.toml`):**
```toml
[input]
targets = "panel.fa"

[output]
output_dir = "results/"
output_format = "tsv,vcf"

[calling]
min_confidence = 0.99
max_vaf = 0.35
```

**Command:**
```bash
kam run \
  --config configs/tumour-informed.toml \
  --r1 new_plasma_R1.fq.gz \
  --r2 new_plasma_R2.fq.gz \
  --target-variants new_patient_mutations.vcf \
  --output-dir results/new_patient/
```

---

### Stage-by-stage pipeline (Nextflow integration)

Run each stage individually. Useful for debugging, benchmarking per-stage performance, or integration with a workflow manager.

```bash
# Stage 1: assemble molecules
kam assemble --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
  --output molecules.bin

# Stage 2: build k-mer index
kam index --input molecules.bin --targets panel.fa \
  --output index.bin --sv-junctions sv_junctions.fa

# Stage 3: walk paths
kam pathfind --index index.bin --targets panel.fa \
  --output paths.bin

# Stage 4: call variants
kam call --paths paths.bin --output variants.tsv \
  --output-format tsv,vcf --max-vaf 0.35
```
