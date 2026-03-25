# kam Pipeline: Overview

## What kam does

kam is an alignment-free variant detection pipeline for duplex UMI sequencing. It takes paired-end
FASTQ files from a UMI panel and calls somatic variants without ever aligning reads to a reference
genome.

The pipeline has four stages run sequentially:

```
FASTQ → [assemble] → molecules → [index] → k-mer index → [pathfind] → scored paths → [call] → variants
```

Each stage is described in its own document:

- `assemble.md` — read parsing, UMI grouping, consensus calling
- `index.md` — k-mer extraction and molecule evidence indexing
- `pathfind.md` — de Bruijn graph construction and path enumeration
- `call.md` — statistical variant calling and filtering
- `output.md` — all output files and their fields

---

## Why alignment-free

Alignment-based callers (Mutect2, Sentieon TNscope) map reads to a reference genome and then
identify positions where reads disagree with the reference. This approach has two limitations
for ctDNA variant calling at low VAF:

1. **Alignment artefacts**: at 0.1–2% VAF, the variant signal is weak. Soft-clipped reads,
   misaligned indels, and mapping quality penalties can all suppress or distort the signal.
2. **Information loss**: standard aligners discard molecule identity. Two reads from the same
   original DNA fragment are treated as independent evidence. Duplex information (two reads from
   opposite strands of the same fragment) is not preserved through the alignment BAM.

kam builds a de Bruijn graph directly from the observed k-mers in each target window. Variant
detection is then a graph path enumeration problem: find all sequences through the graph that
differ from the reference. This approach:

- Never discards read pairs based on mapping quality
- Preserves molecule identity (which reads came from the same fragment)
- Preserves duplex status (which molecules have coverage on both strands)
- Counts evidence in molecules, not reads

---

## The molecule as the atomic unit

Every counting and filtering decision in kam operates at the molecule level, not the read level.

A **molecule** in kam is a single DNA fragment as identified by its UMI pair. It may have:
- One or more reads on the forward strand (simplex forward)
- One or more reads on the reverse strand (simplex reverse)
- Reads on both strands (duplex)

The consensus sequence is computed per strand and then merged if both strands are present.
Duplex consensus is the highest quality evidence: sequencing errors on one strand are corrected
by the other strand.

When kam reports `NALT=3`, it means 3 molecules support the alternate allele, not 3 reads.
At 2M read pairs with 2% VAF, this typically means 3 original DNA fragments carried the variant.

---

## Chemistry support

### Twist UMI duplex (default)

The Twist UMI duplex panel uses the `5M2S+T` read structure on both R1 and R2:

- Positions 0–4: UMI (5 random bases)
- Positions 5–6: skip/spacer (2 monotemplate bases, used for QC)
- Positions 7+: template (the actual genomic sequence)

Each read pair carries two UMIs: the R1 UMI (from the forward strand) and the R2 UMI (from the
reverse strand). The **canonical UMI pair** is the lexicographically smaller of `umi_r1 + umi_r2`
and `umi_r2 + umi_r1`. This makes the identifier strand-agnostic: both the forward and reverse
reads from the same fragment produce the same canonical UMI.

With 5-base UMIs (4^5 = 1,024 possible values per strand), the canonical pair space has at most
1,024 × 1,025 / 2 = 524,800 unique entries. At 2M reads, approximately 67% of this space is
occupied, giving roughly 350K unique molecules.

### Configurable chemistry

The UMI length and skip length are configurable via `config.toml` or CLI flags. This supports
other protocols:

- **12 bp UMI, no skip** (simplex protocols): `umi_length = 12`, `skip_length = 0`, `duplex = false`
- **Custom duplex**: any UMI length with `duplex = true`

See the `examples/` directory for ready-to-use config files and `configuration.md` for the
full reference.

---

## Configuration

All pipeline settings can be specified in a TOML config file:

```toml
[input]
r1 = "sample_R1.fq.gz"
r2 = "sample_R2.fq.gz"
targets = "panel_targets.fa"

[output]
output_dir = "results"
```

Run the pipeline with:

```bash
kam run --config config.toml
```

CLI flags always take precedence over config file values. See `configuration.md` for the full
reference of all sections and fields.

---

## Variant types

kam classifies each detected variant into one of 11 types:

| Type | Size | Description |
|------|------|-------------|
| `SNV` | 1 bp | Single nucleotide variant |
| `MNV` | 2+ bp | Multi-nucleotide variant (same length, multiple positions differ) |
| `Insertion` | 1–49 bp | Short insertion |
| `Deletion` | 1–49 bp | Short deletion |
| `Complex` | — | Complex rearrangement |
| `LargeDeletion` | ≥50 bp | Large deletion, detected via SV junction sequences |
| `TandemDuplication` | ≥50 bp | Tandem duplication, detected via SV junction sequences |
| `Inversion` | ≥50 bp | Inversion, detected via SV junction sequences |
| `Fusion` | — | Junction between two distinct genomic loci (translocation or gene fusion) |
| `InvDel` | ≥50 bp | Inversion with flanking deletion |
| `NovelInsertion` | ≥50 bp | Novel insertion without sequence similarity to flanking reference |

SNVs, MNVs, insertions, and deletions are detected from the standard graph walk. SVs require
either SV junction sequences (`sv_junctions` in config) or fusion target sequences
(`fusion_targets` in config).

---

## Targeted vs discovery mode

kam has two operating modes:

**Discovery mode** (default): finds all variant paths through each target window that pass
statistical quality filters. Reports everything passing filters as `PASS`. Produces 35–72 false
positive calls per sample at 2M reads due to background biological signal in cfDNA (germline
variants, clonal haematopoiesis, somatic mutations in normal tissue).

**Monitoring mode** (`target_variants` in config, or `--target-variants` at CLI): restricts
output to calls that exactly match a pre-specified list of `(CHROM, POS, REF, ALT)` tuples.
All other quality-passing calls are relabelled `NotTargeted`. Produces zero false positives for
a correctly specified truth set. This mode is appropriate for serial ctDNA monitoring where the
somatic variant panel is known from a prior tissue biopsy.

---

## Key numbers at 2M reads, 2% VAF

| Metric | Value |
|--------|-------|
| Molecules assembled | ~300–430K |
| Duplex fraction | 15–21% |
| Targets | 375 |
| Targets with alt paths found | ~270–280 |
| Overall sensitivity | 61–66% |
| SNV sensitivity | 69–80% |
| Indel sensitivity | 32–39% |
| False positives (discovery mode) | 35–72 |
| False positives (monitoring mode) | 0 |
| Runtime (single core) | 19–35s |
| Peak RSS | 1.5–2.0 GB |

---

## Sensitivity ceiling at 2M reads

The sensitivity limit at 2M reads is not molecule dropout. UMI saturation means that additional
reads beyond 2M give diminishing returns: 10M reads produces only 18% more unique molecules than
2M reads (467K vs 351K). The two structural limits are:

1. **Indel end-anchor displacement**: for an indel of length ℓ, the variant path through the
   de Bruijn graph is ℓ k-mers longer or shorter than the reference path. Reads must span to the
   shifted anchor position. At 2% VAF, 61–68% of indels are missed for this reason.
2. **max_paths=100 walk budget**: the DFS walk stops after finding 100 paths. For targets with
   complex graphs (high-complexity sequence, repeat regions), the budget fills with spurious paths
   before the true variant is found.

See `docs/research/depth_saturation_investigation.md` and
`docs/research/target_length_investigation.md` for the full analysis.
