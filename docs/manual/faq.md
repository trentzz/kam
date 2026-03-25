# FAQ

---

## What chemistries are supported?

kam supports any chemistry where both reads in a pair carry a UMI at a fixed position at the
start of the read. The UMI length, skip (spacer) length, and duplex mode are all configurable
via `config.toml` or CLI flags.

Two ready-to-use config presets are in the `examples/` directory:

- `twist-umi-duplex.toml`: Twist Biosciences duplex UMI panel. 5 bp UMI, 2 bp skip, duplex
  (both strands sequenced).
- `simplex-umi-12bp.toml`: a simplex protocol with a 12 bp UMI and no skip region.

For other protocols, set `umi_length`, `skip_length`, and `duplex` in the `[chemistry]` section
to match your read structure. The template starts at position `umi_length + skip_length`.

Protocols without UMIs are not supported. The UMI is required for molecule identity and duplex
assembly.

---

## What variant types can kam detect?

kam detects 11 variant types:

| Type | Minimum requirements |
|------|---------------------|
| `SNV` | Standard graph walk |
| `MNV` | Standard graph walk |
| `Insertion` (< 50 bp) | Standard graph walk |
| `Deletion` (< 50 bp) | Standard graph walk |
| `Complex` | Standard graph walk |
| `LargeDeletion` (≥ 50 bp) | `sv_junctions` file with deletion junction sequences |
| `TandemDuplication` (≥ 50 bp) | `sv_junctions` file with duplication junction sequences |
| `Inversion` (≥ 50 bp) | `sv_junctions` file with inversion junction sequences |
| `InvDel` (≥ 50 bp) | `sv_junctions` file with inversion-deletion junction sequences |
| `NovelInsertion` (≥ 50 bp) | Standard graph walk with extended path length |
| `Fusion` | `fusion_targets` file with synthetic breakpoint junction sequences |

SNVs, MNVs, and short indels are detected from the standard de Bruijn graph walk and do not
require any additional input files. Large SVs require pre-computed junction sequences.

---

## How does tumour-informed monitoring mode work?

In tumour-informed monitoring mode, you provide a VCF of known somatic variants (typically from
a matched tissue biopsy or prior diagnostic test). kam calls variants normally, then applies a
post-calling filter:

- Calls that exactly match a `(CHROM, POS, REF, ALT)` entry in the VCF are marked `PASS`.
- All other quality-passing calls are relabelled `NotTargeted` and excluded from the PASS set.

This produces near-zero false positives. The background biological signal in cfDNA (germline
variants, clonal haematopoiesis) that produces 35–72 PASS calls per sample in discovery mode is
completely suppressed.

To use monitoring mode, set `target_variants` in the `[input]` section of your config:

```toml
[input]
target_variants = "tumour_biopsy.vcf"
```

Or pass `--target-variants tumour_biopsy.vcf` at the CLI.

The `ti_position_tolerance` setting (default: 0) allows a position offset in base pairs for
matching. For large SVs with breakpoint ambiguity, set this to a small value such as 5.

---

## What is the minimum VAF kam can detect?

The practical minimum VAF depends on read depth and the statistical thresholds.

At the default thresholds (`min_confidence = 0.99`, `min_alt_molecules = 2`):

- At 2M reads, Twist UMI duplex, ~350K unique molecules: the minimum detectable VAF is
  approximately 0.3–0.5% (requiring 2 alt molecules out of ~500–700 ref + alt in the target
  window).
- At 5M reads (~750K unique molecules): the minimum detectable VAF drops to approximately
  0.1–0.2%.

In monitoring mode with a known variant, a single duplex-confirmed molecule (n_alt = 1,
n_duplex_alt ≥ 1) can be reported as PASS. This enables detection at even lower VAF when the
variant is pre-specified.

The key constraint is the number of unique molecules per target, not the number of reads. With
5 bp UMIs and 1,024 possible values per strand, the UMI space is small. At 2M reads, it is ~67%
saturated. Adding more reads yields diminishing returns in new unique molecules.

---

## How does kam compare to alignment-based tools?

**Alignment-based callers** (Mutect2, Sentieon TNscope, GATK HaplotypeCaller):
- Map reads to a reference genome.
- Call variants at positions where reads disagree with the reference.
- Discard molecule identity: two reads from the same original fragment are treated as
  independent evidence.
- Cannot natively use duplex information from UMI sequencing.

**kam**:
- Never aligns reads to a reference.
- Builds a de Bruijn graph from k-mers in each target window.
- Counts evidence in molecules (original DNA fragments), not reads.
- Preserves duplex status: knows which molecules have coverage on both strands.
- Produces molecule-level VAF estimates with calibrated credible intervals.

The practical difference: kam's sensitivity is limited by the number of unique molecules (UMI
saturation), while alignment-based callers are limited by mapping quality and alignment artefacts.
For liquid biopsy panel sequencing with 5 bp UMIs and 2M reads, kam achieves 61–66% overall
sensitivity at 2% VAF in discovery mode. Alignment-based callers typically achieve 70–85% at the
same VAF, primarily because they handle indels better. At 0.1% VAF, kam in monitoring mode
achieves near-zero false positives, which is difficult for alignment-based callers without
tumour-informed filtering.

---

## Can I use kam without UMIs?

No. The UMI is required for molecule identity. The pipeline assembles molecules by grouping read
pairs with the same canonical UMI, clustering by Hamming distance, and splitting by endpoint
fingerprint. Without a UMI, there is no way to identify which reads came from the same original
DNA fragment, and duplex assembly is impossible.

For non-UMI sequencing, alignment-based callers with PCR duplicate marking are the appropriate
approach.

---

## Can I use kam for germline variant calling?

kam is designed for somatic variant calling at low VAF. The statistical model uses a background
error rate of 1e-4 and a uniform Beta prior, which is appropriate for rare somatic variants.

For germline variant calling (expected VAF 0.5 or 1.0), the model will assign very high
confidence to any variant with adequate molecule support, and the `max_vaf` filter can be used
to exclude somatic calls. However, the pipeline does not produce genotype information (GT field)
or follow VCF conventions for germline variant representation.

If your goal is germline variant calling from UMI data, consider using alignment-based tools
with UMI-aware duplicate marking (e.g. fgbio, UMI-tools) followed by standard germline callers.

---

## Why does the output include filtered variants?

By default, kam writes all variant calls to the output file, including those that fail a filter.
The `filter` column shows the reason. This allows downstream analysis to apply different
thresholds without rerunning the pipeline.

To extract only PASS calls from TSV output:

```bash
awk 'NR==1 || $NF=="PASS"' variants.tsv > pass_variants.tsv
```

In VCF output, the FILTER column contains `PASS` for passing calls and the filter label
(e.g. `StrandBias`) for filtered calls. Standard VCF tools respect this convention.

---

## What does `n_duplex_alt = 0` mean for a PASS call?

`n_duplex_alt` is the number of duplex molecules whose k-mers overlap the variant-specific (alt)
k-mers. A value of 0 means no duplex molecule has reads on both strands at the variant site.

This is common at low VAF with 2M reads. At 2M reads and 2% VAF with Twist UMI duplex
chemistry, approximately 25–36% of genuine variants have zero variant-specific duplex. The call
is still valid: it is supported by simplex molecules from both strands (reflected in `n_simplex_alt`
and the strand bias p-value).

A non-zero `n_duplex_alt` provides stronger evidence, but is not required by default. The
`min_alt_duplex` threshold (default: 0) controls whether duplex confirmation is required for PASS.

---

## How do I detect structural variants?

SV detection requires additional input files:

1. **SV junction sequences** (`sv_junctions`): a FASTA of synthetic sequences that span the
   breakpoints of LargeDeletion, TandemDuplication, Inversion, and InvDel events. These are
   typically pre-computed for a panel and are specific to the structural variants being monitored.

2. **Fusion target sequences** (`fusion_targets`): a FASTA of synthetic sequences encoding
   breakpoint junctions between two genomic loci. Each entry ID must follow the format
   `{name}__{chromA}:{startA}-{endA}__{chromB}:{startB}-{endB}__fusion`.

Set these in the `[input]` section of your config. The SV-specific calling thresholds
(`sv_min_confidence`, `sv_min_alt_molecules`, `sv_strand_bias_threshold`) can be tuned
independently from the SNV/indel thresholds.

If you do not have pre-computed SV junction sequences, SVs will not be detected. The standard
de Bruijn graph walk does not find structural variants that span beyond the target window.
