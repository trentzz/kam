# Interpreting kam Results

## Overview

kam produces variant calls with full statistical evidence. This document explains what each
output field means, how to read the QC JSON files, what patterns to expect, and what to watch
for.

---

## TSV/VCF columns explained

### target_id

The identifier of the target window where the variant was found. This comes directly from the
FASTA header of the target sequence file.

When target IDs follow the coordinate format `chrN:START-END`, VCF output will contain proper
genomic coordinates. Otherwise, VCF output uses `target_id` as CHROM and position 1.

Examples:
- `chr17:7674220-7674320` — a 100 bp window on chr17.
- `TP53_exon7` — a named target; VCF will use this as CHROM.

### variant_type

All 11 variant type values:

| Value | Meaning |
|-------|---------|
| `SNV` | Single nucleotide variant (same length, 1 base differs) |
| `MNV` | Multi-nucleotide variant (same length, 2+ bases differ) |
| `Insertion` | Short insertion (<50 bp, alt longer than ref) |
| `Deletion` | Short deletion (<50 bp, alt shorter than ref) |
| `Complex` | Complex rearrangement (not currently produced by the walker) |
| `LargeDeletion` | Deletion ≥ 50 bp |
| `TandemDuplication` | Insertion ≥ 50 bp, alt is a tandem repeat of nearby reference |
| `Inversion` | Same-length rearrangement, alt is reverse complement of ref segment |
| `NovelInsertion` | Insertion ≥ 50 bp, not a tandem repeat |
| `InvDel` | Inversion combined with flanking deletion |
| `Fusion` | Junction between two distinct genomic loci |

### vaf, vaf_ci_low, vaf_ci_high

`vaf` is the point estimate of variant allele frequency: `n_alt / (n_ref + n_alt)`.

`vaf_ci_low` and `vaf_ci_high` are the 2.5th and 97.5th percentiles of the Beta posterior
`Beta(n_alt + 1, n_ref + 1)`. This is a 95% Bayesian credible interval.

**How to read confidence intervals**: wide intervals mean few molecules. At `n_alt=2`,
`n_ref=500` (VAF=0.4%), the 95% CI is approximately [0.05%, 1.4%]. At `n_alt=10`, `n_ref=500`
(VAF=2%), the CI is approximately [1.0%, 3.5%].

A call with `vaf_ci_low` close to zero is uncertain. Even if it passes the confidence filter,
treat it as provisional.

### n_molecules_ref, n_molecules_alt

Molecule counts, not read counts. Each molecule is one original DNA fragment as identified
by its canonical UMI pair. Two reads from the same fragment count as one molecule.

`n_molecules_ref` is the minimum molecule count across k-mers in the reference path. For
targets with many reads but a narrow bottleneck in the graph, this can be lower than expected.

`n_molecules_alt` for SNVs and short indels is the minimum molecule count across k-mers in
the alt path. For SV types it is the mean over variant-specific k-mers only (see
`guides/structural-variants.md`).

### n_duplex_alt, n_simplex_alt

`n_duplex_alt` is the count of alt-supporting molecules that are duplex: the same original
fragment was sequenced on both strands, and both strands independently support the alt allele.

`n_simplex_alt` is the count of alt-supporting molecules on a single strand only
(`n_molecules_alt - n_duplex_alt`).

Duplex support is strong evidence of a real variant. Both strands of the same DNA fragment
must independently contain the variant. The probability of a concordant sequencing error on
both strands is ~10⁻⁶.

At 2M reads with 15–21% duplex fraction, approximately 25–36% of genuine 2% VAF variants
have zero variant-specific duplex coverage. Absence of duplex support does not rule out a
real variant, but its presence strongly confirms one.

### strand_bias_p

Fisher's exact test two-tailed p-value for strand bias. The test checks whether the alt allele
is distributed evenly between forward and reverse strands, relative to the reference allele.

A value of 1.0 means perfect balance. Values below 0.01 (the default threshold) cause the
call to be labelled `StrandBias`.

Genuine somatic variants appear on both strands because the variant is present in the original
double-stranded DNA. Single-strand artefacts (oxidative damage, end-repair errors) are
strand-specific.

**Exception**: inversions are structurally strand-biased due to directional path walking in
the de Bruijn graph. The `sv_strand_bias_threshold` defaults to 1.0 (disabled) for SV types.

### confidence

The posterior probability that the observed alt molecules are real signal rather than
sequencing error. Computed as a logistic function of the binomial log-likelihood ratio between
the observed VAF and the background error rate (default 1e-4).

Values close to 1.0 indicate that the observed molecule count is far above what background
error alone would produce. At `n_alt=10`, `n_total=500` (VAF=2%), `confidence ≈ 1.0`.
At `n_alt=2`, `n_total=500` (VAF=0.4%), `confidence` may drop toward threshold.

The default minimum threshold is 0.99 for SNVs/indels and 0.95 for SV types.

### filter

Seven possible values:

| Value | Meaning |
|-------|---------|
| `PASS` | All quality filters passed. The call is reported as a genuine variant. |
| `LowConfidence` | Posterior confidence below threshold, or fewer alt molecules than `min_alt_molecules`. |
| `StrandBias` | Alt allele is predominantly on one strand (Fisher p < `strand_bias_threshold`). |
| `LowDuplex` | Variant-specific duplex count below `min_alt_duplex` (disabled by default). |
| `HighVaf` | VAF exceeds `max_vaf` (disabled by default). Used to exclude germline heterozygous variants. |
| `CollisionRisk` | UMI collision risk too high. Not currently assigned by the pipeline. |
| `NotTargeted` | Tumour-informed monitoring mode: the call does not match `--target-variants`. |

Filters are applied in this order: molecule count, duplex count, strand bias, confidence,
max VAF, tumour-informed. The first failing filter determines the label.

All calls are included in the TSV output regardless of filter. To get only PASS calls:

```bash
awk '$14 == "PASS"' variants.tsv
```

---

## How to read the QC JSON files

### assembly_qc.json

Written by `kam assemble`. Key fields:

| Field | What it tells you |
|-------|-------------------|
| `n_read_pairs_processed` | Total read pairs in the FASTQ |
| `n_read_pairs_passed` | Pairs that passed QC (good UMI quality, minimum length) |
| `n_molecules` | Unique DNA fragments (UMI pairs) assembled |
| `n_duplex` | Molecules with coverage on both strands |
| `n_simplex_fwd` / `n_simplex_rev` | Single-strand molecules by orientation |
| `n_umi_collisions_detected` | Fragments sharing a UMI pair (collision-aware resolution applied) |

**What to check**: `n_molecules` should be 200K–500K at 2M reads. A value below 100K suggests
a library quality issue. `n_duplex / n_molecules` is the duplex fraction; typical values are
15–21% at 2M reads.

### index_qc.json

Written by `kam index`. Key fields:

| Field | What it tells you |
|-------|-------------------|
| `n_targets` | Number of target sequences |
| `n_kmers_observed` | Total on-target k-mers found in the data |
| `n_molecules_indexed` | Molecules with at least one on-target k-mer |

**What to check**: `n_molecules_indexed` should be close to `n_molecules` from assembly QC.
A large gap indicates poor target coverage in the library.

### pathfind_qc.json

Written by `kam pathfind`. Key fields:

| Field | What it tells you |
|-------|-------------------|
| `n_targets` | Targets queried |
| `no_anchor` | Targets with no anchor k-mers found (not enough molecule coverage) |
| `no_paths` | Targets where no path was found (graph walk failed) |
| `ref_only` | Targets with a reference path but no alt paths |
| `with_variants` | Targets with at least one non-reference path found |

**What to check**: `no_anchor` should be low (< 5% of targets). A high `no_anchor` count means
molecules are not reaching these target windows, which is a library quality or target design
issue. `ref_only` is normal: most targets have no somatic variant.

### call_qc.json

Written by `kam call`. Key fields:

| Field | What it tells you |
|-------|-------------------|
| `n_variants_called` | Total alt paths passed to the caller |
| `n_pass` | Calls that received `PASS` |
| `n_filtered` | Calls with any non-PASS filter |

**What to check**: `n_pass` in discovery mode is typically 35–72 per 2M-read sample. Values
above 100 may indicate a high-noise library. Values of 0 may indicate a configuration error.
In monitoring mode, `n_pass` should equal the number of target variants that are present
in the sample.

---

## Common patterns

### What a true positive looks like

```
target_id       variant_type    vaf       vaf_ci_low  vaf_ci_high  n_ref  n_alt  n_duplex  n_simplex  sb_p    conf      filter
chr17:7674220   SNV             0.019800  0.011200    0.032500     495    10     3         7          0.842   0.999999  PASS
```

Characteristics:
- `filter = PASS`
- `n_duplex_alt >= 1`
- `strand_bias_p > 0.1` (no strand bias)
- `confidence >= 0.999`
- `vaf_ci_low` clearly above zero
- Both `n_ref` and `n_alt` are plausible given the VAF

### What a false positive looks like

```
target_id       variant_type    vaf       vaf_ci_low  vaf_ci_high  n_ref  n_alt  n_duplex  n_simplex  sb_p    conf      filter
chr7:117548800  SNV             0.003900  0.000500    0.014000     507    2      0         2          0.031   0.991     StrandBias
```

Characteristics:
- Single-strand alt support (`n_duplex_alt = 0`)
- `strand_bias_p = 0.031`, below threshold
- Wide credible interval
- `filter = StrandBias`

### What a missed variant looks like

When a variant is missed, there is no entry in the output. Common reasons:

- **Low molecule count at the target**: check `no_anchor` in pathfind QC.
- **Variant at the target edge**: indels near the edge shift the anchor k-mer outside the window.
- **Very low VAF**: at < 0.1% VAF, stochastic molecule dropout removes the alt path.

A reference-only entry (all filters `LowConfidence`, `n_alt = 0` or 1) at the correct
`target_id` indicates the target was queried but no alt path was found or passed filters.

### Interpreting VAF confidence intervals

The 95% credible interval reflects molecule count uncertainty. Rules of thumb:

| n_alt | Approximate CI width |
|-------|----------------------|
| 2 | ~0 to ~3× the point VAF |
| 5 | ~0.4× to ~2× the point VAF |
| 10 | ~0.6× to ~1.7× the point VAF |
| 20 | ~0.7× to ~1.4× the point VAF |
| 50+ | CI narrows to ±30% of point VAF |

Use the CI when deciding whether two samples have the same VAF. If the intervals overlap, the
VAF difference is not statistically significant at 95% confidence.

### When to trust discovery mode vs tumour-informed mode

**Use discovery mode** when:
- You are profiling a new sample without prior tissue biopsy.
- You want to find all variants, including unexpected ones.
- You accept 35–72 background calls per sample and will filter downstream.

**Use tumour-informed mode** when:
- You have a matched tissue biopsy VCF.
- You are doing serial monitoring (tracking tumour burden over time).
- You need near-zero false positives for clinical reporting.

Discovery mode's 35–72 background calls per sample come from germline variants not excluded
by `--max-vaf`, clonal haematopoiesis, and somatic mosaicism in normal tissue. These are
real biological signals, not sequencing artefacts. Running a 0% VAF (no tumour DNA) control
on the same panel establishes the expected background for your specific panel design.

---

## Red flags

### Very low on-target molecule count

If `n_molecules_indexed` in `index_qc.json` is less than 50K, the library has poor target
capture. Most targets will have too few molecules to call variants reliably. Check:
- Library preparation quality (ligation efficiency, capture efficiency).
- Target FASTA correctness (headers, sequence orientation).
- Whether reads are being dropped by UMI quality filters (check `n_low_umi_quality` in
  `assembly_qc.json`).

### High PASS rate in 0% VAF controls

If a 0% VAF control sample (no tumour DNA spiked in) shows more than 100 PASS calls in
discovery mode, investigate:
- **Contamination**: cross-sample contamination from a high-tumour-burden sample.
- **Background somatic signal**: clonal haematopoiesis in the normal donor produces real
  somatic variants that pass all filters.
- **Artefact**: systematic library prep artefact generating false alt paths. Check whether
  the same `target_id` recurs across all samples.

Use the 0% control PASS list as a blocklist for your panel.

### Strand bias on all variants

If most PASS calls have `strand_bias_p < 0.05`, a systematic strand-specific artefact is
present. Common causes:
- Oxidative damage (G→T transversions predominantly on one strand).
- End-repair artefact introducing false insertions at read ends.
- Chemistry mismatch (wrong `--chemistry` setting).

Lower `--strand-bias-threshold` to 0.001 to filter these more aggressively, at the cost of
some sensitivity for real low-VAF variants.
