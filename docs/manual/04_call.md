# kam Pipeline: Stage 4 — Call

## What this stage does

The call stage converts scored paths into variant calls. For each target, it takes the reference path and every non-reference path and computes:

- A VAF (variant allele frequency) estimate with 95% credible interval.
- A posterior probability that the variant is real.
- A strand bias test.
- A filter label.

All arithmetic uses molecule counts, not read counts.

---

## Input

The input is a `ScoredPathRecord` file from the pathfind stage. Records are grouped by `target_id`. Within each target, one path is marked `is_reference = true`; all others are candidate variant paths.

For each alt path, `call_variant` is invoked with the reference path evidence and the alt path evidence.

---

## VAF estimation

### Point estimate

The VAF point estimate is:

```
VAF = n_alt / (n_ref + n_alt)
```

where `n_alt = alt_evidence.min_molecules` and `n_ref = ref_evidence.min_molecules`.

The minimum molecule count is used rather than a raw sum because the minimum represents the weakest supported k-mer in the path. A path with one very low-support k-mer should be treated as having low support overall, not averaged away.

### 95% Beta credible interval

The Beta posterior for VAF uses a conjugate prior Beta(1, 1) (uniform) updated with the molecule counts:

```
posterior = Beta(k + 1, M - k + 1)
```

where `k = n_alt` and `M = n_ref + n_alt`.

The 95% credible interval is the 2.5th and 97.5th percentiles of this Beta distribution:
- `vaf_ci_low = Beta(k+1, M-k+1).quantile(0.025)`
- `vaf_ci_high = Beta(k+1, M-k+1).quantile(0.975)`

At k=5, M=1000 (VAF=0.5%): the CI is approximately [0.16%, 1.16%], reflecting high uncertainty at very low counts. At k=20, M=1000 (VAF=2%): the CI is approximately [1.2%, 3.0%].

---

## Strand bias test

The strand bias test checks whether the alt-supporting molecules are distributed unevenly between forward and reverse strands relative to the reference molecules.

The test uses Fisher's exact test on a 2×2 contingency table:

```
           fwd strand   rev strand
alt            a              b
ref            c              d
```

where:
- `a = alt_evidence.min_simplex_fwd`
- `b = alt_evidence.min_simplex_rev`
- `c = ref_evidence.min_simplex_fwd`
- `d = ref_evidence.min_simplex_rev`

The two-tailed p-value is computed by summing the hypergeometric probabilities of all 2×2 tables with fixed marginals that are at least as extreme as the observed table (i.e., whose hypergeometric probability is ≤ the observed probability).

A small p-value indicates that the alt allele is predominantly on one strand. For genuine somatic variants, the alt should appear on both strands. All-one-strand calls are almost always sequencing artefacts.

The default threshold is `strand_bias_threshold = 0.01`. A call with `strand_bias_p < 0.01` is labelled `StrandBias`.

---

## Posterior confidence

The posterior confidence is a Bayesian estimate of the probability that the observed alt molecules are real (signal) rather than background sequencing error.

The computation uses a log-likelihood ratio between two hypotheses:

- **Signal**: VAF = k/M (the observed rate is the true rate).
- **Background**: VAF = `background_error_rate` (the observed rate is background error).

```
log_LR = log P(k | VAF=k/M) - log P(k | VAF=background_error_rate)
```

where `P(k | VAF)` is the binomial likelihood of observing `k` successes in `M` trials with probability `VAF`.

With equal priors, the posterior is:

```
confidence = exp(log_LR) / (1 + exp(log_LR))
```

The default `background_error_rate = 1e-4` (1 in 10,000). This is typical of Phred 40 sequencing error after consensus calling.

At k=10, M=1000 (VAF=1%): `confidence ≈ 1.0` because 10 molecules is far more than expected from background error.
At k=2, M=1000 (VAF=0.2%): `confidence` may drop below 0.99 depending on the exact background rate.

A call with `confidence < min_confidence` is labelled `LowConfidence`.

---

## Variant classification

Variants are classified by comparing the lengths and bases of the reference and alt sequences:

| Classification | Condition |
|---------------|-----------|
| `SNV` | Same length, exactly 1 position differs |
| `MNV` | Same length, >1 positions differ |
| `Insertion` | Alt is longer than ref |
| `Deletion` | Alt is shorter than ref |
| `Complex` | Complex rearrangement (not currently produced by the walker) |

The sequences compared are the full reconstructed path sequences, not left-normalised VCF alleles. Left normalisation to VCF coordinates is applied separately when writing VCF output (see `05_outputs.md`).

---

## Filter assignment

Filters are applied in a fixed priority order. The first failing filter determines the label.

### 1. Molecule count check

`n_alt < min_alt_molecules`

The default `min_alt_molecules = 2` requires at least 2 molecules supporting the alt allele. A single-molecule call (n_alt = 1) is very likely a sequencing error unless duplex-confirmed.

**Exception**: a single-molecule call is accepted if all of the following are true:
1. `n_alt == 1`.
2. `min_alt_duplex_for_single` is set (default: `Some(1)`).
3. The variant-specific duplex count `n_duplex_alt >= min_alt_duplex_for_single`.

Duplex confirmation on a single molecule is strong evidence. If both strands of the same fragment independently show the variant, the probability of a concordant error on both strands is the product of both per-strand error rates (e.g., ~10⁻⁶). This exception enables detection at very low depth (e.g., low-coverage targets) without raising the overall false positive rate.

Calls failing this check: `LowConfidence`.

### 2. Variant-specific duplex check

`n_duplex_alt < min_alt_duplex`

The default `min_alt_duplex = 0` disables this check. When set to 1, every PASS call must have at least one duplex molecule whose k-mers overlap the variant site.

This filter is intentionally disabled by default. At 2M reads (6–7% duplex fraction at the variant site), approximately 25–36% of genuine 2% VAF variants have zero variant-specific duplex. Enabling this filter would remove too many true positives at this depth. Consider enabling it at higher duplex fractions (≥15%) or when running at higher depth (≥5M reads).

Calls failing this check: `LowDuplex`.

### 3. Strand bias check

`strand_bias_p < strand_bias_threshold`

The default `strand_bias_threshold = 0.01`. Fisher's exact test p-value below this threshold flags the call as strand-biased.

Calls failing this check: `StrandBias`.

### 4. Confidence check

`confidence < min_confidence`

The default `min_confidence = 0.99`. Calls with posterior probability below this threshold are too ambiguous.

Calls failing this check: `LowConfidence`.

### 5. Maximum VAF check

`vaf > max_vaf`

Disabled by default (`max_vaf = None`). When set, calls with VAF above this threshold are labelled `HighVaf`. This is intended to exclude germline heterozygous variants (VAF ≈ 0.5) from a somatic calling run. A typical ctDNA setting is `--max-vaf 0.35`.

Calls failing this check: `HighVaf`.

### 6. PASS

All filters passed. The call is a genuine somatic variant.

---

## Filter labels

| Label | Meaning |
|-------|---------|
| `PASS` | All filters passed |
| `LowConfidence` | Posterior confidence below threshold, or molecule count too low |
| `StrandBias` | Strand bias Fisher p-value below threshold |
| `LowDuplex` | Variant-specific duplex count below threshold |
| `HighVaf` | VAF exceeds maximum threshold (likely germline) |
| `CollisionRisk` | UMI collision risk too high (not currently assigned by the pipeline) |
| `NotTargeted` | Monitoring mode: call does not match the target variants list |

---

## Tumour-informed tumour-informed mode

When `--target-variants VCF` is provided, a post-calling filter is applied:

1. Load the truth VCF as a set of `(CHROM, POS, REF, ALT)` tuples.
2. For each PASS call, check whether it matches an entry in the truth set.
3. If it does not match: relabel as `NotTargeted`.

Only calls that exactly match a pre-specified truth variant remain `PASS`. All other quality-passing calls are suppressed. This produces zero false positives for a correctly specified truth set (e.g., from a matched tissue biopsy).

The tumour-informed filter is the reason kam achieves near-zero false positives in tracking mode. The background biology in cfDNA (germline variants, clonal haematopoiesis, somatic mosaicism in normal tissue) produces ~35–72 PASS calls per sample in discovery mode. Monitoring mode filters all of them.

---

## CallerConfig summary

| Field | Default | CLI flag |
|-------|---------|---------|
| `min_confidence` | 0.99 | `--min-confidence` |
| `strand_bias_threshold` | 0.01 | `--strand-bias-threshold` |
| `min_alt_molecules` | 2 | `--min-alt-molecules` |
| `min_alt_duplex` | 0 | `--min-alt-duplex` |
| `min_alt_duplex_for_single` | Some(1) | Not exposed |
| `max_vaf` | None | `--max-vaf` |
| `background_error_rate` | 1e-4 | Not exposed |

---

## Call QC

The call QC JSON (`call_qc.json`) records:

| Field | Meaning |
|-------|---------|
| `n_variants_called` | Total variant paths processed |
| `n_pass` | Calls with `PASS` label |
| `n_filtered` | Calls with any non-PASS label |

---

## Output

Variants are written in one or more formats (see `05_outputs.md`). The default is TSV. Multiple formats can be requested with a comma-separated list: `--output-format tsv,vcf`.
