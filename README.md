# kam

Alignment-free variant detection for Twist duplex UMI sequencing.

kam replaces the four-tool chain HUMID + Jellyfish + km + kmtools with a single
Rust binary that preserves molecule-level evidence from raw FASTQ reads through
to variant calls. It supports two operating modes: somatic discovery and
tumour-informed monitoring.

---

## Features

- **Alignment-free**: no reference alignment required at call time.
- **Molecule-aware**: every k-mer carries a `MoleculeEvidence` record with
  distinct molecule count, duplex support, and strand breakdown.
- **Two modes**: discovery (all alt paths) and tumour-informed monitoring
  (`--target-variants`).
- **Tumour-informed monitoring**: precision 1.0 at all VAF levels. Background
  cfDNA variants are eliminated because they do not match the expected somatic
  alleles.
- **Fast**: 16--22 seconds per sample on a single core at 2 M read pairs.

---

## Benchmark results

Evaluated on the Twist cfDNA Pan-Cancer Reference Standard v2 (24 samples,
375 truth variants, 3 concentrations, 8 VAF levels).
Configuration: `--max-vaf 0.35 --min-family-size 2 --target-variants`.

| Concentration | VAF | Sensitivity | SNV sens | Indel sens | Precision |
|---------------|-----|-------------|----------|-----------|-----------|
| 15 ng         | 2%  | 61.3%       | 80.0%    | 38.8%     | 1.0       |
| 30 ng         | 2%  | 59.2%       | 77.1%    | 37.6%     | 1.0       |
| 5 ng          | 2%  | 51.7%       | 68.8%    | 31.2%     | 1.0       |
| 15 ng         | 0.5%| 40.0%       | 52.7%    | 24.7%     | 1.0       |
| 30 ng         | 0.5%| 46.1%       | 61.0%    | 28.2%     | 1.0       |
| 0% VAF (all)  | —   | 0 FPs       | —        | —         | —         |

Runtime: 16--22 s per sample. Peak RSS: 1.8--2.0 GB.

---

## Installation

### From crates.io

```sh
cargo install kam-bio
```

This installs the `kam` binary.

### From source

Requires Rust 1.75 or later.

```sh
cargo build --release
```

The binary is at `target/release/kam`.

---

## Usage

### End-to-end run (recommended)

```sh
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets targets_100bp.fa \
  --output-dir results/
```

#### Tumour-informed monitoring mode

Supply a VCF of expected somatic variants. Only calls matching an entry in the
target set are reported; all other calls are labelled `NotTargeted`.

```sh
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets targets_100bp.fa \
  --output-dir results/ \
  --max-vaf 0.35 \
  --min-family-size 2 \
  --target-variants known_variants.vcf
```

### Individual pipeline stages

The pipeline can also be run stage by stage:

```sh
# 1. Assemble molecules from raw reads.
kam assemble --r1 R1.fastq.gz --r2 R2.fastq.gz --output molecules.bin

# 2. Build a k-mer index against target sequences.
kam index --molecules molecules.bin --targets targets.fa --output index.bin

# 3. Walk de Bruijn graph paths.
kam pathfind --index index.bin --targets targets.fa --output paths.bin

# 4. Call variants from scored paths.
kam call --paths paths.bin --targets targets.fa --output calls.vcf
```

### Key options

| Flag | Default | Description |
|------|---------|-------------|
| `--min-family-size N` | 1 | Minimum reads per UMI family. Set to 2 to remove singletons. |
| `--max-vaf F` | — | Discard calls above this VAF (removes germline heterozygotes). |
| `--min-alt-molecules N` | 2 | Minimum alt molecules to emit a call. |
| `--min-confidence F` | 0.99 | Minimum posterior confidence. |
| `--target-variants VCF` | — | Enable tumour-informed monitoring mode. |

---

## Chemistry

Designed for Twist UMI duplex chemistry (`5M2S+T` read structure on both R1
and R2). The 5 bp random UMI, 2 bp skip, and template are extracted and used
to group read families and identify duplex pairs.

---

## Architecture

```
kam-core      — shared types: Molecule, ConsensusRead, MoleculeEvidence
kam-assemble  — molecule assembly from raw FASTQ (replaces HUMID)
kam-index     — k-mer indexing with molecule provenance (replaces Jellyfish)
kam-pathfind  — de Bruijn graph construction and path walking (replaces km)
kam-call      — statistical variant calling and tumour-informed filtering
kam           — CLI binary wiring all stages together
```

---

## Development

Run tests and quality checks before committing:

```sh
cargo test
cargo clippy -- -D warnings
cargo fmt -- --check
```

---

## Licence

MIT
