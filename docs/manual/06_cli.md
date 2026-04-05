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

For most uses, `kam run` is the correct entry point. The individual subcommands are exposed for debugging, benchmarking individual stages, and integration with Nextflow.

---

## kam run

Run the complete pipeline: assemble → index → pathfind → call.

```
kam run --r1 <R1.fastq.gz> --r2 <R2.fastq.gz> --targets <targets.fa> --output-dir <DIR> [OPTIONS]
```

### Required arguments

| Flag | Type | Description |
|------|------|-------------|
| `--r1` | Path | R1 input FASTQ file (gzip-compressed or plain) |
| `--r2` | Path | R2 input FASTQ file (gzip-compressed or plain) |
| `--targets` | Path | Target sequences FASTA file |
| `--output-dir` | Path | Directory for all pipeline outputs |

### Assembly options

| Flag | Default | Description |
|------|---------|-------------|
| `--chemistry` | `twist-umi-duplex` | Chemistry preset. Currently only `twist-umi-duplex` is supported. Determines the UMI length (5 bp), skip length (2 bp), and template start position (7 bp). |
| `--min-umi-quality` | 20 | Minimum Phred quality for UMI bases. Read pairs with any UMI base below this quality are dropped before grouping. Set to 0 to disable. |
| `--min-family-size` | 1 | Minimum reads per UMI family. Families with fewer reads are discarded. Setting to 2 eliminates true singletons but reduces sensitivity at low depth. |
| `--min-template-length` | (none) | Minimum template length in bases. Templates shorter than this on either R1 or R2 are dropped. Useful for filtering very short inserts. Not set by default; the k-mer extraction stage handles minimum effective length implicitly. |

### Indexing options

| Flag | Default | Description |
|------|---------|-------------|
| `-k` / `--kmer-size` | 31 | K-mer length. Must be 1–31. Larger k improves specificity (fewer spurious paths) but requires reads that fully span the k-mer window. At k=31, reads shorter than 31 bp produce no k-mers. Benchmark data suggests k=31 is the sweet spot for 100 bp targets at 2M reads. |

### Calling options

| Flag | Default | Description |
|------|---------|-------------|
| `--min-confidence` | 0.99 | Minimum posterior probability for a PASS call. The posterior is derived from a binomial likelihood ratio between the observed VAF and the background error rate. Values below this threshold produce LowConfidence calls. Lower values increase sensitivity but increase false positives. |
| `--strand-bias-threshold` | 0.01 | Fisher's exact test p-value cutoff for strand bias. Calls with p below this threshold are labelled StrandBias. Genuine somatic variants appear on both strands. Artefacts (e.g. oxidative damage, end-repair artefacts) are strand-specific. The default 0.01 is a standard threshold; lowering it reduces false positives from strand-specific artefacts at the cost of some sensitivity. |
| `--min-alt-molecules` | 2 | Minimum molecules supporting the alt allele. A single-molecule call (n_alt=1) is accepted only if duplex-confirmed (see `min_alt_duplex_for_single`). Setting to 3 or higher increases specificity but reduces sensitivity for low-VAF variants. |
| `--min-alt-duplex` | 0 | Minimum variant-specific duplex molecules. Counts only duplex molecules whose k-mers are at the variant site (alt-specific k-mers), not at the flanking anchor k-mers. Default 0 disables this filter. Set to 1 to require at least one duplex confirmation at the variant site. Recommended only at duplex fractions ≥15% or depth ≥5M reads. |
| `--max-vaf` | (none) | Maximum VAF for a PASS call. Calls with VAF above this are labelled HighVaf. For ctDNA somatic calling, germline heterozygous variants have VAF ≈ 0.5. Setting `--max-vaf 0.35` eliminates them while retaining all low-VAF somatic calls. |
| `--target-variants` | (none) | VCF file of expected somatic variants. When set, only calls matching a `(CHROM, POS, REF, ALT)` entry in this VCF are marked PASS. All other quality-passing calls are labelled NotTargeted. This is tumour-informed monitoring mode. Use it when the somatic variant panel is known from a prior tissue biopsy. Produces near-zero false positives. |

### Output options

| Flag | Default | Description |
|------|---------|-------------|
| `--output-format` | `tsv` | Comma-separated list of output formats. Supported values: `tsv`, `csv`, `json`, `vcf`. When multiple formats are given, each is written to `<output_base>.<ext>` using the output file base name. |
| `--qc-output` | (none) | Path for a merged QC JSON file. Not currently produced; per-stage QC JSONs are written to the output directory. |

### Logging options

| Flag | Default | Description |
|------|---------|-------------|
| `--log-dir` | (none) | Directory for structured log output. Not currently implemented. |
| `--log` | (none, repeatable) | Enable specific log sinks. Can be repeated: `--log umi --log family`. Not currently implemented. |

### Performance options

| Flag | Default | Description |
|------|---------|-------------|
| `--threads` | (none) | Number of worker threads. Not currently implemented; the pipeline runs single-threaded. |

---

## kam assemble

Parse paired FASTQ files and produce a molecules bincode file.

```
kam assemble --r1 <R1.fastq.gz> --r2 <R2.fastq.gz> --output <molecules.bin> [OPTIONS]
```

### Required arguments

| Flag | Description |
|------|-------------|
| `--r1` | R1 input FASTQ |
| `--r2` | R2 input FASTQ |
| `--output` | Output file path for assembled molecules |

### Optional arguments

Same assembly options as `kam run` above:

| Flag | Default |
|------|---------|
| `--chemistry` | `twist-umi-duplex` |
| `--min-umi-quality` | 20 |
| `--min-family-size` | 1 |
| `--min-template-length` | (none) |
| `--log-dir` | (none) |
| `--log` | (none, repeatable) |
| `--threads` | (none) |

The assemble stage writes `assemble_qc.json` to the same directory as `--output`.

---

## kam index

Build a k-mer index from assembled molecules.

```
kam index --input <molecules.bin> --targets <targets.fa> --output <index.bin> [OPTIONS]
```

### Required arguments

| Flag | Description |
|------|-------------|
| `--input` | Molecules file (output of `kam assemble`) |
| `--targets` | Target sequences FASTA |
| `--output` | Output k-mer index file |

### Optional arguments

| Flag | Default | Description |
|------|---------|-------------|
| `-k` / `--kmer-size` | 31 | K-mer length (1–31) |

The index stage writes `index_qc.json` to the same directory as `--output`.

---

## kam pathfind

Walk de Bruijn graph paths through the k-mer index for each target.

```
kam pathfind --index <index.bin> --targets <targets.fa> --output <paths.bin>
```

### Required arguments

| Flag | Description |
|------|-------------|
| `--index` | K-mer index file (output of `kam index`) |
| `--targets` | Target sequences FASTA |
| `--output` | Output scored paths file |

There are no optional arguments currently exposed at the CLI for pathfind. The walk configuration (`max_paths = 100`, `max_path_length = 150`) is compiled in.

The pathfind stage writes `pathfind_qc.json` to the same directory as `--output`.

---

## kam call

Call variants from scored paths.

```
kam call --paths <paths.bin> --output <variants.tsv> [OPTIONS]
```

### Required arguments

| Flag | Description |
|------|-------------|
| `--paths` | Scored paths file (output of `kam pathfind`) |
| `--output` | Output variant calls file |

### Optional arguments

Same calling options as `kam run` above:

| Flag | Default |
|------|---------|
| `--output-format` | `tsv` |
| `--min-confidence` | (CallerConfig default: 0.99) |
| `--strand-bias-threshold` | (CallerConfig default: 0.01) |
| `--min-alt-molecules` | (CallerConfig default: 2) |
| `--min-alt-duplex` | (CallerConfig default: 0) |
| `--max-vaf` | (none) |
| `--target-variants` | (none) |

The call stage writes `call_qc.json` to the same directory as `--output`.

---

## Output format detail

### Multiple formats

Passing a comma-separated list to `--output-format` writes each format to a separate file. The output path base (without extension) is used:

```
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

## Examples

### Full pipeline, discovery mode

```bash
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --output-dir ./results \
  --output-format vcf
```

### Full pipeline, tumour-informed monitoring mode

```bash
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --output-dir ./results \
  --target-variants tumour_biopsy.vcf \
  --output-format vcf \
  --max-vaf 0.35
```

### Run with duplex confirmation required

```bash
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --output-dir ./results \
  --min-alt-duplex 1 \
  --output-format tsv,vcf
```

### Stage-by-stage pipeline

```bash
kam assemble \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --output molecules.bin

kam index \
  --input molecules.bin \
  --targets panel_targets.fa \
  --output index.bin \
  -k 31

kam pathfind \
  --index index.bin \
  --targets panel_targets.fa \
  --output paths.bin

kam call \
  --paths paths.bin \
  --output variants.vcf \
  --output-format vcf \
  --min-confidence 0.99 \
  --max-vaf 0.35 \
  --target-variants tumour_biopsy.vcf
```

---

## Notes

### File formats

All intermediate binary files (`*.bin`) are `bincode`-serialised with a header identifying the file type and the kam version. Passing a file produced by one version of kam to a different version will fail with an error if the serialisation format changed.

### Determinism

All pipeline stages produce deterministic output for a given input, regardless of the number of threads. Internal HashMap ordering is not relied upon for output: molecules are sorted by UMI hash, variants are sorted by target ID before writing.

### Error messages

All errors print to stderr. The exit code is 1 on any error, 0 on success. Stage progress is reported to stderr with `[stage]` prefixes, e.g. `[assemble] molecules=351341 duplex=66552`.
