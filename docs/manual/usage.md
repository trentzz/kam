# kam Usage Guide

## Getting started

### Install

```bash
cargo install --path .
```

Or from the workspace root:

```bash
cargo build --release
# Binary is at: target/release/kam
```

### First run

The minimum required inputs are two FASTQ files and a target FASTA:

```bash
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --output-dir results/
```

This runs the full pipeline (assemble → index → pathfind → call) and writes `variants.tsv` to
`results/`. Per-stage QC JSON files are also written to `results/`.

### Using a config file vs CLI flags

Most settings can be provided in a TOML config file with `--config`. This is useful for
reproducible runs and Nextflow integration.

```toml
# kam.toml
[input]
r1 = "sample_R1.fastq.gz"
r2 = "sample_R2.fastq.gz"
targets = "panel_targets.fa"

[output]
dir = "results/"
format = ["tsv", "vcf"]

[calling]
max_vaf = 0.35
min_confidence = 0.99
```

```bash
kam run --config kam.toml
```

CLI flags override config file values when both are present.

---

## Common workflows

### Discovery mode: finding all variants

Discovery mode (the default) reports all variant paths that pass quality filters. No prior
knowledge of the variants is required.

```bash
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --output-dir results/ \
  --max-vaf 0.35 \
  --output-format tsv,vcf
```

`--max-vaf 0.35` excludes germline heterozygous variants (VAF ≈ 0.5) from the output. This
is recommended for all somatic ctDNA runs.

In discovery mode, expect 35–72 PASS calls per sample at 2M reads. These include germline
variants not filtered by `--max-vaf`, clonal haematopoiesis calls, and somatic mosaic variants
in normal tissue. Run a 0% VAF control sample on the same panel to establish the biological
background.

### Monitoring mode: tracking known variants

Tumour-informed monitoring restricts output to variants that exactly match a pre-specified
VCF. All other quality-passing calls are labelled `NotTargeted`.

```bash
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --output-dir results/ \
  --max-vaf 0.35 \
  --target-variants biopsy_somatic.vcf \
  --output-format tsv,vcf
```

This produces near-zero false positives. Use it for serial ctDNA monitoring after a matched
tissue biopsy.

The `biopsy_somatic.vcf` must contain (CHROM, POS, REF, ALT) tuples matching the expected
somatic variants. CHROM should match the `target_id` format in your target FASTA (either
`chrN:START-END` or a named target ID).

### SV detection: adding junction k-mers

Large deletions, tandem duplications, and novel insertions are detected automatically — no extra
flags needed. Junction sequences are only required for **inversions and InvDel events**, where
the breakpoint k-mers are absent from the standard panel targets. See
`guides/structural-variants.md` for how to create junction sequences.

```bash
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --sv-junctions sv_junctions.fa \
  --output-dir results/ \
  --max-vaf 0.35 \
  --output-format tsv,vcf
```

SVs appear in the VCF with symbolic allele notation (`<DEL>`, `<DUP>`, `<INV>`, `<INS>`,
`<INVDEL>`) and additional INFO tags `SVTYPE` and `SVLEN`.

### Fusion detection: adding fusion targets

Fusion detection requires synthetic junction sequences via `--fusion-targets`. See
`guides/fusions.md` for the FASTA header format and how to create fusion targets.

```bash
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --fusion-targets fusion_targets.fa \
  --output-dir results/ \
  --output-format tsv,vcf
```

Fusions are written as paired BND records in VCF output.

### Everything at once

Run SNV/indel, SV, and fusion detection in a single pass:

```bash
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --sv-junctions sv_junctions.fa \
  --fusion-targets fusion_targets.fa \
  --target-variants biopsy_somatic.vcf \
  --output-dir results/ \
  --max-vaf 0.35 \
  --output-format tsv,vcf
```

---

## Input requirements

### FASTQ format

- Paired-end reads (two files: R1 and R2).
- Gzip-compressed (`.fastq.gz`) or plain (`.fastq`).
- Read structure: Twist UMI duplex chemistry (`5M2S+T` on both R1 and R2). The first 5 bases
  are the UMI, bases 6–7 are the skip spacer, and base 8 onward is the template.

### Target FASTA

The target FASTA contains one entry per panel target. Each entry is the reference sequence
for that target window.

**Header format options:**

Use genomic coordinate headers for proper VCF output:
```
>chr17:7674220-7674320
ACGT...
```

Or named headers (VCF will use the name as CHROM, POS=1):
```
>TP53_exon7
ACGT...
```

**Recommended target design:**

- Window length: 100–150 bp. Shorter windows miss indels near the edges; longer windows
  increase graph complexity and runtime.
- Centre the region of interest: place the expected variant position in the middle 40 bp
  of the window.
- For indels: ensure the window is wide enough that the shifted anchor k-mer (offset by the
  indel length) is still within the window.

### SV junction FASTA

One entry per SV. Each entry is the breakpoint-spanning sequence. See
`guides/structural-variants.md` for construction details.

No special header format is required. The header is used as the `target_id` in the output.

### Fusion target FASTA

One entry per fusion. The header must follow the exact format:

```
>{name}__{chromA}:{startA}-{endA}__{chromB}:{startB}-{endB}__fusion
```

See `guides/fusions.md` for details.

### Target variants VCF (tumour-informed mode)

A standard VCF with (CHROM, POS, REF, ALT) columns. CHROM must match the `target_id` in
your target FASTA. Only the four core VCF columns are used; other fields are ignored.

---

## Output files

| File | Contents |
|------|----------|
| `variants.tsv` | All variant calls (PASS and filtered). Tab-separated, one call per line. |
| `variants.vcf` | Same calls in VCF 4.3 format (if `--output-format` includes `vcf`). |
| `assembly_qc.json` | Molecule assembly statistics. |
| `index_qc.json` | K-mer indexing statistics. |
| `pathfind_qc.json` | Path finding statistics (anchor hits, paths found). |
| `call_qc.json` | Variant calling statistics (total, PASS, filtered). |

The key result file is `variants.tsv` (or `variants.vcf`). All other files are QC diagnostics.

To get only PASS calls from the TSV:

```bash
awk 'NR==1 || $14 == "PASS"' variants.tsv > pass_calls.tsv
```

To count PASS calls:

```bash
awk '$14 == "PASS"' variants.tsv | wc -l
```

---

## Tips

### Start with defaults, then tune

Default settings are calibrated for 2M read pairs, 100 bp targets, and 2% VAF somatic calls.
Start with defaults and adjust only if:
- Sensitivity is too low (try lowering `--min-confidence` to 0.95 or `--min-alt-molecules` to 1).
- False positive rate is too high (try raising `--min-alt-duplex` to 1, or enable
  `--target-variants` monitoring mode).

### Always run 0% VAF control samples

Run the pipeline on a control sample with no tumour DNA. This establishes the panel-specific
background: germline variants not filtered by `--max-vaf`, clonal haematopoiesis, and
systematic artefacts. Use the PASS call list from the control as a blocklist for your study
samples.

### Use tumour-informed mode for clinical monitoring

If you have a matched tissue biopsy VCF, use `--target-variants`. Discovery mode produces
35–72 background calls per sample; monitoring mode produces zero false positives for a
correctly specified truth set. For ctDNA monitoring applications, the monitoring mode result
is the one to report.

### Check QC JSON for library quality

Before trusting variant calls, check:
- `n_molecules` in `assembly_qc.json`: should be 200K–500K at 2M reads.
- `n_duplex / n_molecules` in `assembly_qc.json`: duplex fraction should be 15–21%.
- `no_anchor` in `pathfind_qc.json`: should be < 5% of targets.

Low molecule counts or high `no_anchor` rates indicate a library quality issue that will limit
sensitivity regardless of caller settings.

### Understand the sensitivity limits

At 2M reads, kam achieves:
- 90% SNV sensitivity at 0.15% VAF, 100% at 0.2% VAF (discovery mode).
- 90% indel sensitivity at 0.2% VAF, 100% at 0.25% VAF (discovery mode).
- 100% SV sensitivity at 0.25% VAF (discovery mode with junction k-mers).

These numbers come from synthetic titration benchmarks (varforge). Real samples may differ
due to non-uniform target coverage, GC bias, and panel design.

The primary limit below 0.2% VAF is molecule dropout: at very low VAF, the expected number
of alt molecules is below 1, and stochastic sampling means some samples have zero alt
molecules even though the variant is present. Increasing sequencing depth beyond 2M reads
gives diminishing returns due to UMI saturation (10M reads adds only 18% more unique molecules).

### Runtime and memory

At 2M reads with 375 targets, a single-core run completes in 19–35 seconds and uses 1.5–2.0 GB
peak RSS. Parallelism is not yet implemented; plan resource allocation accordingly for Nextflow.
