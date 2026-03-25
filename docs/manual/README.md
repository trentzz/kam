# kam User Manual

kam is an alignment-free variant detection pipeline for duplex UMI sequencing. It takes
paired-end FASTQ files from a UMI panel and calls somatic variants without aligning reads to a
reference genome.

Target use case: detecting somatic variants at 0.1–2% VAF in liquid biopsy ctDNA panel samples,
with or without a matched tumour sample for filtering.

---

## Quick Start

### Install

```bash
cargo install kam-bio
```

Or build from source:

```bash
git clone https://github.com/trentzz/kam
cd kam
cargo build --release
```

### Minimal example

```bash
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --output-dir results/
```

This runs the full pipeline and writes `results/variants.tsv`.

### Using a config file

```bash
kam run --config config.toml
```

Config files keep all settings in one place and are recommended for reproducible runs. See
`examples/twist-umi-duplex.toml` for the Twist UMI duplex preset and
`examples/simplex-umi-12bp.toml` for a simplex 12 bp UMI protocol.

---

## Table of Contents

### Getting Started

- [Configuration](configuration.md) — full config.toml reference
- [CLI reference](modules/overview.md#targeted-vs-discovery-mode) — quick orientation

### Pipeline Modules

Each stage of the pipeline is described in detail:

- [Overview](modules/overview.md) — pipeline architecture, variant types, discovery vs monitoring
- [Assemble](modules/assemble.md) — read parsing, UMI grouping, consensus calling
- [Index](modules/index.md) — k-mer extraction, molecule evidence index, SV junction k-mers
- [Pathfind](modules/pathfind.md) — de Bruijn graph, DFS walk, SV fallback paths
- [Call](modules/call.md) — VAF estimation, strand bias, filters, fusion calling
- [Output](modules/output.md) — TSV/CSV/JSON/VCF formats, QC JSON schemas, fusion BND

### Reference

- [Troubleshooting](troubleshooting.md) — common problems and solutions
- [FAQ](faq.md) — frequently asked questions
