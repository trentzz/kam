# Discovery Mode Precision Study

**Date:** 2026-04-13
**Benchmark:** 07-snvindel-ml-boost-v1 (24-sample titration)
**Question:** How close can discovery mode get to tumour-informed (TI) mode precision, and does the ML model help?

---

## Background

TI mode achieves precision = 1.000 by construction: only calls whose (CHROM, POS, REF, ALT) exactly matches a truth VCF are passed. Discovery mode applies no such restriction and calls all variants that pass statistical filters. This study characterises the precision gap, tests whether the `twist-duplex-v1` ML model can close it, and identifies the best practical operating point for discovery-mode calling.

---

## Conditions

Five conditions were run across all 24 titration samples (3 ng × 8 VAF), 2M reads per sample:

| Label | Flags | Output |
|-------|-------|--------|
| A — Discovery | _(none)_ | `titration_2mreads_disc.tsv` |
| B — Discovery + ML | `--ml-model twist-duplex-v1` | `titration_2mreads_disc_ml.tsv` |
| C — Discovery + ML + germline | `--ml-model twist-duplex-v1 --max-vaf 0.35` | `titration_2mreads_disc_ml_maxvaf35.tsv` |
| D — Discovery + ML + germline + duplex | `--ml-model twist-duplex-v1 --max-vaf 0.35 --min-alt-duplex 1` | `titration_2mreads_disc_ml_maxvaf35_d1.tsv` |
| E — Discovery + germline only | `--max-vaf 0.35` | `titration_2mreads_disc_maxvaf35.tsv` |

Compared against the existing TI-only run (`titration_2mreads_ti.tsv`).

Per-sample variants TSVs (with `ml_prob`) saved to `results/tsvs_disc_ml/` for condition B.

---

## Key Finding: The ML Model Does Not Work in Discovery Mode

Every PASS call in discovery mode — both true positives and false positives — receives `ml_prob` in [0.9, 1.0]. The threshold cannot be set to separate the two populations because they are indistinguishable to the model.

| | TP calls | FP calls |
|---|---|---|
| ml_prob in [0.9, 1.0] | 1974 / 1974 (100%) | 1590 / 1590 (100%) |
| ml_prob >= 0.5 | 1974 / 1974 (100%) | 1590 / 1590 (100%) |

The precision-recall curve is flat at every threshold from 0.0 to 0.99: raising the ML threshold removes nothing.

Conditions A and B are identical in every metric: same FP count, same sensitivity, same precision at every sample and VAF level. Adding `--ml-model twist-duplex-v1` in discovery mode has zero effect.

**Why this happens.** The statistical caller requires posterior confidence ≥ 0.99 by default. Any call that passes this threshold already exhibits the profile that the ML model learned to associate with true positives: adequate molecule support, low strand bias, and high confidence. The model was trained on synthetic data where the only PASS calls were real variants, so it has no reference distribution for real background errors at high confidence. In simulation, sequencing errors are filtered before the calling stage; in real data they are not. The model cannot distinguish a background C→T deamination that happens to pass the stats filter from a genuine somatic SNV.

Conditions B and C also give identical FP counts: all ML-filtered FP counts in condition C equal the germline-filtered FP counts in condition E. This confirms that the ML filter removes zero calls on top of what `--max-vaf 0.35` already removes.

---

## Sensitivity is Identical in Discovery and TI Mode

TI mode does not improve sensitivity. It only restricts which PASS calls are counted as PASS, eliminating FPs but leaving TPs unaffected. The underlying statistical calling is identical between discovery and TI mode.

| ng | VAF | Discovery sens | TI sens |
|----|-----|---------------|---------|
| 15ng | 0.5% | 0.424 | 0.424 |
| 15ng | 1.0% | 0.576 | 0.576 |
| 15ng | 2.0% | 0.627 | 0.627 |
| 30ng | 2.0% | 0.608 | 0.608 |
| 5ng  | 2.0% | 0.539 | 0.539 |

Sensitivity differences only appear when `--min-alt-duplex 1` is added (condition D), which does affect PASS/FAIL decisions for low-duplex calls.

---

## Precision and FP Counts by Condition

### 15 ng

| VAF% | Disc | +Germ | +ML+Germ | +ML+G+Dup | TI | Disc FP | +Germ FP | +ML+G+D FP | TI FP |
|------|------|-------|----------|-----------|-----|---------|----------|------------|-------|
| 0.001 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 70 | 61 | 8 | 0 |
| 0.01  | 0.015 | 0.017 | 0.017 | 0.000 | 1.000 | 67 | 59 | 20 | 0 |
| 0.1   | 0.256 | 0.282 | 0.282 | 0.529 | 1.000 | 64 | 56 | 8 | 0 |
| 0.25  | 0.591 | 0.618 | 0.618 | 0.836 | 1.000 | 76 | 68 | 12 | 0 |
| 0.5   | 0.685 | 0.713 | 0.713 | 0.933 | 1.000 | 73 | 64 | 6 | 0 |
| 1.0   | 0.755 | 0.776 | 0.776 | 0.924 | 1.000 | 70 | 62 | 12 | 0 |
| 2.0   | 0.693 | 0.709 | 0.709 | 0.883 | 1.000 | 104 | 96 | 28 | 0 |
| 0% (neg) | 0.013 | 0.015 | 0.015 | 0.000 | 1.000 | 74 | 66 | 7 | 0 |

### 30 ng

| VAF% | Disc | +Germ | +ML+Germ | +ML+G+Dup | TI | Disc FP | +Germ FP | +ML+G+D FP | TI FP |
|------|------|-------|----------|-----------|-----|---------|----------|------------|-------|
| 0.001 | 0.024 | 0.027 | 0.027 | 0.000 | 1.000 | 81 | 73 | 13 | 0 |
| 0.01  | 0.024 | 0.027 | 0.027 | 0.000 | 1.000 | 81 | 73 | 16 | 0 |
| 0.1   | 0.287 | 0.310 | 0.310 | 0.421 | 1.000 | 77 | 69 | 11 | 0 |
| 0.25  | 0.578 | 0.605 | 0.605 | 0.909 | 1.000 | 78 | 70 | 3 | 0 |
| 0.5   | 0.686 | 0.706 | 0.706 | 0.855 | 1.000 | 85 | 77 | 11 | 0 |
| 1.0   | 0.720 | 0.738 | 0.738 | 0.936 | 1.000 | 86 | 78 | 7 | 0 |
| 2.0   | 0.752 | 0.772 | 0.772 | 0.943 | 1.000 | 75 | 67 | 9 | 0 |
| 0% (neg) | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 86 | 78 | 6 | 0 |

### 5 ng

| VAF% | Disc | +Germ | +ML+Germ | +ML+G+Dup | TI | Disc FP | +Germ FP | +ML+G+D FP | TI FP |
|------|------|-------|----------|-----------|-----|---------|----------|------------|-------|
| 0.001 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 34 | 26 | 4 | 0 |
| 0.01  | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 42 | 34 | 10 | 0 |
| 0.1   | 0.054 | 0.062 | 0.062 | 0.118 | 1.000 | 53 | 45 | 15 | 0 |
| 0.25  | 0.384 | 0.431 | 0.431 | 0.478 | 1.000 | 45 | 37 | 12 | 0 |
| 0.5   | 0.648 | 0.700 | 0.700 | 0.870 | 1.000 | 38 | 30 | 7 | 0 |
| 1.0   | 0.761 | 0.794 | 0.794 | 0.922 | 1.000 | 47 | 39 | 9 | 0 |
| 2.0   | 0.838 | 0.866 | 0.866 | 0.956 | 1.000 | 39 | 31 | 8 | 0 |
| 0% (neg) | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 43 | 34 | 7 | 0 |

---

## FP Breakdown: What Each Filter Addresses

At 15 ng with no spiked-in variants (0% VAF), the 74 PASS calls are all false positives. Their composition can be estimated from the condition comparisons:

- **Germline variants** (removed by `--max-vaf 0.35`): ~8 FPs per sample. These are heterozygous germline variants at VAF ≈ 0.5. The filter removes them without affecting sensitivity because no truth variant has VAF ≥ 0.35.
- **Low-duplex background errors** (removed by `--min-alt-duplex 1`): ~57–62 FPs per sample. Background sequencing errors typically produce one or two simplex alt molecules with no duplex confirmation. This filter removes them but also removes real low-VAF calls with NDUPALT=0 (approximately 25–49% of true positives at 2M reads, depending on depth and VAF).
- **Residual FPs** (persist through all filters in condition D): 3–28 per sample. These have NDUPALT ≥ 1, VAF ≤ 0.35, and high confidence. They are likely recurrent error hotspots at specific k-mer contexts, or low-level PCR artefacts that generate enough reads to create a duplex pair by chance.

The ML model addresses none of these categories because it assigns probability ≥ 0.9 to all of them.

---

## Sensitivity Cost of Each Filter

Sensitivity is unchanged between conditions A, B, C, and E — all discovery runs without `--min-alt-duplex 1` have identical TP counts. Only condition D (with `--min-alt-duplex 1`) costs sensitivity, and that cost was characterised in `duplex-filter-analysis.md`.

Summary of sensitivity vs precision trade-off at 2% VAF:

| Condition | 15ng prec | 15ng sens | 15ng FP | 30ng prec | 30ng sens | 30ng FP |
|-----------|-----------|-----------|---------|-----------|-----------|---------|
| Discovery | 0.693 | 0.627 | 104 | 0.752 | 0.608 | 75 |
| + Germline filter | 0.709 | 0.627 | 96 | 0.772 | 0.608 | 67 |
| + Germline + Duplex | 0.883 | 0.565 | 28 | 0.943 | 0.395 | 9 |
| TI mode | 1.000 | 0.627 | 0 | 1.000 | 0.608 | 0 |

TI mode achieves maximum precision at zero sensitivity cost. No discovery configuration matches this.

---

## Operating Point Recommendations

**For somatic monitoring (known variants, ctDNA):** Use TI mode (`--target-variants`). It is the only configuration that achieves precision = 1.000 with no sensitivity loss. Discovery mode is not appropriate for monitoring applications.

**For discovery calling (unknown variants, panel sequencing):** Always include `--max-vaf 0.35` to suppress germline heterozygous variants. This removes ~8 FPs per sample at no sensitivity cost. Do not add `--ml-model` — it provides no benefit at 2M reads.

**If high precision is required in discovery mode:** Add `--min-alt-duplex 1` on top of `--max-vaf 0.35`. This achieves precision of 0.88–0.96 at VAF ≥ 1% but costs 10–35% sensitivity depending on DNA input. The 30 ng input is hit hardest. Only appropriate at sequencing depths ≥ 10M reads where duplex capture is more complete.

**Do not use `--ml-model` in discovery mode.** The `twist-duplex-v1` model provides zero improvement over unfiltered discovery calling. It was trained on synthetic data where all statistical PASS calls are true positives. Real background errors that pass the stats filter are indistinguishable from real variants in the model's feature space. Retraining on real data (true positives from the titration dataset, false positives from the negative controls) is required before the ML model can be useful in discovery mode.

---

## Implication for ML Model Development

The failure of `twist-duplex-v1` in discovery mode points to a clear training data problem rather than a feature problem. The 33 features available (VAF, molecule counts, duplex counts, strand bias, confidence, sequence context) are in principle sufficient to distinguish many background errors from somatic variants: background errors typically have lower duplex support, occur at specific k-mer contexts associated with oxidative or deamination damage, and cluster at known error hotspots.

The model cannot exploit these signals because it was never trained on examples of real background errors. Every training example labelled "negative" in the synthetic dataset is a variant that was not spiked in — it does not resemble a real sequencing artefact. The model has learned a boundary between "called with high confidence" and "not called", not between "somatic" and "artefact".

To make ML useful in discovery mode, the training set needs real FP examples. The titration negative control samples (0% VAF) provide exactly this: ~40–86 PASS calls per sample that are confirmed artefacts (no spiked-in variants). Combined with confirmed TPs from the high-VAF samples, this gives a real-data training set. The estimated total across 3 negative controls is ~200–250 FP examples per ng condition, small but representative of the actual error distribution in this chemistry and panel.
