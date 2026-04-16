# SNV and Indel Detection

## What SNVs and indels are in kam

A **single nucleotide variant (SNV)** is a single-base substitution: the alternate sequence has
the same length as the reference and differs at exactly one position.

An **indel** is an insertion or deletion shorter than 50 bp. Insertions have an alt sequence
longer than the reference; deletions have an alt sequence shorter than the reference.

A **multi-nucleotide variant (MNV)** is a same-length substitution at two or more adjacent
positions. kam detects MNVs but they can be noisy (see limitations below).

kam classifies variants after path reconstruction, comparing the full reference path sequence
against the full alt path sequence. Left-normalisation to VCF coordinates is applied separately
when writing VCF output.

---

## How kam detects SNVs and indels

Detection is graph-based, not alignment-based. The steps are:

1. **Index**: k-mers from each molecule are stored with their molecule provenance (UMI pair,
   duplex status, strand orientation).

2. **Graph walk**: for each target window, the k-mer index is treated as a de Bruijn graph.
   A depth-first search walks all paths from the left anchor to the right anchor. The reference
   path is identified by matching the target reference sequence. Every other path is a candidate
   variant.

3. **Bubble detection**: an SNV or indel creates a bubble in the graph. The reference path and
   the alt path diverge at the variant site and reconnect at a shared anchor k-mer. Short indels
   shift the path by the indel length, so the reconnection anchor is offset relative to the
   reference anchor.

4. **Scoring**: each path is scored by the minimum molecule count across all its k-mers. For
   SNVs, this is the minimum of the per-k-mer molecule counts on the variant-containing k-mers.

5. **Calling**: the caller estimates VAF, computes a posterior confidence, tests for strand bias,
   and applies filters (see `04_call.md` for the full model).

---

## Recommended settings for SNV/indel detection

Default settings work well for most ctDNA panel samples. The table below shows the defaults and
when to change them.

| Setting | Default | When to change |
|---------|---------|---------------|
| `--min-confidence` | 0.99 | Lower to 0.95 for lower-depth targets (reduces specificity) |
| `--min-alt-molecules` | 2 | Raise to 3 for high-noise backgrounds |
| `--min-alt-duplex` | 0 | Set to 1 at duplex fraction ≥15% or depth ≥5M reads |
| `--strand-bias-threshold` | 0.01 | Do not lower; increasing increases false positives |
| `--max-vaf` | (none) | Set to 0.35 to exclude germline heterozygous variants |
| `--target-variants` | (none) | Provide a tumour-biopsy VCF for monitoring mode |

For a standard somatic ctDNA panel run:

```
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --output-dir results/ \
  --max-vaf 0.35 \
  --output-format tsv,vcf
```

For tumour-informed monitoring (near-zero FPs):

```
kam run \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --targets panel_targets.fa \
  --output-dir results/ \
  --max-vaf 0.35 \
  --target-variants biopsy_somatic.vcf \
  --output-format tsv,vcf
```

---

## Config file equivalent

```toml
[input]
r1 = "sample_R1.fastq.gz"
r2 = "sample_R2.fastq.gz"
targets = "panel_targets.fa"
target_variants = "biopsy_somatic.vcf"

[output]
output_dir = "results/"
output_format = "tsv,vcf"

[calling]
max_vaf = 0.35
```

---

## What kam is good at

**High sensitivity from 0.25% VAF**: in benchmarks on synthetic titrations, SNV sensitivity
reaches 90% at 0.15% VAF and 100% at 0.2% VAF in discovery mode. Indel sensitivity reaches
90% at 0.2% VAF and 100% at 0.25% VAF.

**Near-zero false positives in tumour-informed mode**: monitoring mode with `--target-variants`
produces zero false positives for a correctly specified truth set. Discovery mode produces
35–72 background calls per sample (germline variants, clonal haematopoiesis) at 2M reads.

**Molecule-level evidence**: all counts are molecule counts, not read counts. NALT=3 means
3 original DNA fragments, not 3 reads.

**Duplex confirmation**: duplex molecules (both strands of the same fragment) provide
high-confidence evidence. A single duplex-confirmed molecule is accepted even when
`min_alt_molecules=2` cannot be met, because concordant error on both strands is improbable
(error rate ~10⁻⁶).

---

## What kam is not good at

**Very low VAF (<0.1%) without tumour-informed mode**: at 0.05% VAF, SNV sensitivity is 60%
in discovery mode. Confidence filters remove many true positives at this depth.

**MNVs**: MNVs can be reported but are more likely to appear as two co-occurring SNVs rather
than a single MNV call, depending on k-mer length and the distance between the substitutions.
Treat MNV calls with caution and inspect the alt sequence manually.

**Short indels near target ends**: indels near the end of the target window shift the required
anchor k-mer outside the window. Sensitivity for indels within 31 bp of either target edge
is reduced.

**Very short targets**: targets shorter than 2×k (default 62 bp) may not support full path
reconstruction. Recommended minimum target length is 100 bp.

---

## Interpreting SNV/indel results

The key fields in the output:

| Field | What it tells you |
|-------|-------------------|
| `vaf` | Point estimate of variant allele frequency |
| `vaf_ci_low` / `vaf_ci_high` | 95% Beta credible interval. Wide intervals mean few molecules |
| `n_molecules_ref` | Molecules supporting the reference allele at this target |
| `n_molecules_alt` | Molecules supporting the variant allele |
| `n_duplex_alt` | Alt molecules that are duplex (both strands confirm the variant) |
| `n_simplex_alt` | Alt molecules on a single strand only |
| `strand_bias_p` | Fisher p-value for strand bias. Values < 0.01 are suspicious |
| `confidence` | Posterior probability the variant is real. Default threshold: 0.99 |
| `filter` | `PASS`, or the reason the call was filtered |

**A strong true positive** has: `filter=PASS`, `n_duplex_alt >= 1`, `strand_bias_p > 0.1`,
`confidence >= 0.999`, and a `vaf_ci_low` that does not include zero.

**A weak but plausible call** has: `PASS`, `n_duplex_alt=0`, `n_molecules_alt=2`,
`confidence` close to threshold. These warrant validation.

**A likely artefact** has: `filter=StrandBias`, or `filter=LowConfidence` with
`n_molecules_alt=1` and `n_duplex_alt=0`.

---

## Tips for optimising detection

- **Use 100 bp target windows**: targets of 100 bp give the best balance of path complexity and
  read spanning. Shorter windows miss indels; longer windows increase graph complexity.

- **Check `n_molecules_ref` per target**: if many targets show very low reference molecule
  counts (<100), library quality or target design may be limiting sensitivity.

- **Inspect `pathfind_qc.json`**: the `no_anchor` count shows how many targets had no
  anchor k-mers found in the data. This is the primary cause of missed variants.

- **At very low VAF, trust duplex**: if `n_duplex_alt >= 1`, a call with `n_molecules_alt=1`
  at 0.1% VAF is more credible than a call with `n_molecules_alt=3` and `n_duplex_alt=0`.

- **Run 0% VAF controls**: the false positive rate in discovery mode is 35–72 per sample.
  Running a control sample establishes the biological background for your panel.
