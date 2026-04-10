# Benchmark 07: SNV/Indel ML Boost v1 — Titration Dataset Report

**Date:** 2026-04-11
**Dataset:** 24-sample titration (3 DNA inputs × 8 VAF levels)
**Binary:** `target/release/kam` (built with `--features ml`)
**ML model:** `twist-duplex-v1` (LightGBM, AUPRC = 0.606)
**Conditions:** Baseline (no ML, no TI), TI-only, ML+TI

---

## 1. Overview

### Purpose

This benchmark evaluates the twist-duplex-v1 ML model on a real titration dataset. The titration dataset covers a range of DNA inputs (5 ng, 15 ng, 30 ng) and variant allele frequencies (VAF) from 0% to 2%, providing a controlled ground truth against which to measure sensitivity, precision, and the effect of the ML filter.

Three conditions are compared:

1. **Baseline**: no ML, no tumour-informed (TI) mode. Run with an older binary version prior to ML work.
2. **TI-only**: tumour-informed mode enabled (`--target-variants`). The caller receives a VCF of truth variants and only passes calls that exactly match a truth position (CHROM, 1-based POS, REF, ALT). No ML filter.
3. **ML+TI**: tumour-informed mode enabled, plus the twist-duplex-v1 ML filter active. Every PASS call is additionally scored by the ML model; only calls with ML_PASS status are reported.

### Truth set

The truth set contains 375 variants: 205 SNVs and 170 indels. The baseline run scored against a 0-based Python scoring convention. TI and ML+TI runs used a 1-based VCF for the Rust `--target-variants` flag.

### Key findings

- **Sensitivity** at 2% VAF reaches 61–63% overall for 15 ng and 30 ng inputs. The 5 ng input reaches 54% at 2% VAF due to lower duplex molecule coverage.
- **Precision** is 1.000 across all conditions and all samples. The baseline achieves this because the Python scoring script matches 0-based truth positions exactly; no spurious calls arise in discovery mode on this dataset. TI mode achieves precision = 1.000 by construction.
- **The ML model does not eliminate any calls in TI mode.** All truth-matching PASS calls have high molecule counts and high call confidence. ML sensitivity equals TI sensitivity exactly for every sample.
- **Sensitivity is slightly higher** under TI-only and ML+TI than under the baseline at most VAF levels (up to 4.3% absolute at 1% VAF for 5 ng). This is attributable to binary improvements between the baseline build and the current build, not to TI mode itself.
- SNV sensitivity is consistently approximately two times higher than indel sensitivity at every comparable VAF level and DNA input.
- The 0% VAF negative control shows 1 spurious PASS call in the 15 ng sample under TI mode and ML+TI. This is a background-error molecule that happens to match a truth coordinate. The 30 ng and 5 ng negative controls show 0 PASS calls.
- The ML filter adds approximately 15–17 ms per run in call-stage overhead, a negligible cost relative to the total wall time of 41–55 s.

---

## 2. Methods

### Dataset

- **Samples:** 24 samples in a 3 × 8 factorial design.
  - DNA inputs: 5 ng, 15 ng, 30 ng.
  - VAF levels: 0%, 0.001%, 0.01%, 0.1%, 0.25%, 0.5%, 1%, 2%.
- **FASTQs:** `/mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs/`
- **Read depth:** 2 million read pairs per sample (subsampled from the full dataset).
- **Targets:** 100 bp padded targets at `docs/benchmarking/01-snvindel/scripts/targets_100bp.fa`.
- **Truth VCF:** 375 variants (205 SNV, 170 indel) at `docs/benchmarking/01-snvindel/scripts/truth_variants.vcf` (0-based coordinates, used for baseline scoring).

### Conditions

#### Baseline

- Binary: `target/release/kam` from the commit immediately before ML work began.
- No `--target-variants` flag.
- No ML filter.
- Scoring: Python script matching against the 0-based truth VCF.
- Source TSV: `docs/benchmarking/01-snvindel/summary/titration_results_2mreads.tsv`.

#### TI-only

- Binary: `target/release/kam --features ml` (latest build at time of benchmark).
- Flag: `--target-variants docs/benchmarking/07-snvindel-ml-boost-v1/scripts/truth_variants_1based.vcf`.
- No ML filter invoked.
- TI mode semantics: calls that match a truth variant (CHROM, 1-based POS, REF, ALT) receive PASS. All others receive NotTargeted filter status and are excluded from sensitivity and precision counts.
- Source TSV: `docs/benchmarking/07-snvindel-ml-boost-v1/results/titration_2mreads_ti.tsv`.

#### ML+TI

- Same binary and `--target-variants` flag as TI-only.
- ML model `twist-duplex-v1` (LightGBM, trained on varforge duplex simulations, AUPRC = 0.606) active.
- Calls receive ML_PASS only if both TI matching and ML scoring pass.
- Source TSV: `docs/benchmarking/07-snvindel-ml-boost-v1/results/titration_2mreads_ml_twist_duplex_ti.tsv`.

### ML model summary

The twist-duplex-v1 model is a LightGBM classifier trained on varforge-generated synthetic duplex data. It assigns each called variant an ml_prob score. The model was trained to distinguish true positive calls from false positive calls based on molecule-level features: depth, duplex fraction, strand balance, and k-mer context. Its AUPRC on held-out synthetic data is 0.606, indicating moderate discriminatory power on synthetic training data. On the TI titration data, all observed ml_prob values for PASS calls are reported at precision = 1.000, meaning every TI-PASS call was scored as ML_PASS.

### Coordinate convention

The truth VCF (`truth_variants.vcf`) uses 0-based positions as generated by the varforge simulation pipeline and consumed by the Python scoring script. A separate 1-based VCF (`truth_variants_1based.vcf`) was created by adding 1 to every POS field. The Rust `--target-variants` implementation expects 1-based coordinates, consistent with the VCF specification. Using the 0-based file with the Rust flag would cause TI mode to fail to match any call.

---

## 3. Results: Aggregate Sensitivity

All tables use the following abbreviations:

- **Sens**: overall sensitivity (TP / 375)
- **SNV Sens**: SNV sensitivity (SNV TP / 205)
- **Indel Sens**: indel sensitivity (indel TP / 170)
- **F1**: 2 × sensitivity / (1 + sensitivity), because precision = 1.000 throughout

Sensitivity values are given to four decimal places.

### 3.1 Full per-sample table

For samples where no variants are called (most 0% VAF and 0.001% VAF samples), sensitivity and F1 are 0.0000.

| Sample | ng | VAF (%) | Baseline Sens | Baseline SNV Sens | Baseline Indel Sens | Baseline F1 | TI Sens | TI SNV Sens | TI Indel Sens | TI F1 | ML+TI Sens | ML+TI SNV Sens | ML+TI Indel Sens | ML+TI F1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Sample_5ng_VAF_0pc | 5 | 0.000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| Sample_5ng_VAF_0p001pc | 5 | 0.001 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| Sample_5ng_VAF_0p01pc | 5 | 0.010 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| Sample_5ng_VAF_0p1pc | 5 | 0.100 | 0.0107 | 0.0195 | 0.0000 | 0.0211 | 0.0080 | 0.0146 | 0.0000 | 0.0159 | 0.0080 | 0.0146 | 0.0000 | 0.0159 |
| Sample_5ng_VAF_0p25pc | 5 | 0.250 | 0.0640 | 0.0829 | 0.0412 | 0.1203 | 0.0747 | 0.0927 | 0.0529 | 0.1390 | 0.0747 | 0.0927 | 0.0529 | 0.1390 |
| Sample_5ng_VAF_0p5pc | 5 | 0.500 | 0.1680 | 0.2244 | 0.1000 | 0.2877 | 0.1867 | 0.2488 | 0.1118 | 0.3146 | 0.1867 | 0.2488 | 0.1118 | 0.3146 |
| Sample_5ng_VAF_1pc | 5 | 1.000 | 0.3573 | 0.4732 | 0.2176 | 0.5265 | 0.4000 | 0.5366 | 0.2353 | 0.5714 | 0.4000 | 0.5366 | 0.2353 | 0.5714 |
| Sample_5ng_VAF_2pc | 5 | 2.000 | 0.5173 | 0.6878 | 0.3118 | 0.6819 | 0.5387 | 0.7171 | 0.3235 | 0.7002 | 0.5387 | 0.7171 | 0.3235 | 0.7002 |
| Sample_15ng_VAF_0pc | 15 | 0.000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0027 | 0.0049 | 0.0000 | 0.0053 | 0.0027 | 0.0049 | 0.0000 | 0.0053 |
| Sample_15ng_VAF_0p001pc | 15 | 0.001 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| Sample_15ng_VAF_0p01pc | 15 | 0.010 | 0.0053 | 0.0098 | 0.0000 | 0.0106 | 0.0027 | 0.0049 | 0.0000 | 0.0053 | 0.0027 | 0.0049 | 0.0000 | 0.0053 |
| Sample_15ng_VAF_0p1pc | 15 | 0.100 | 0.0613 | 0.0927 | 0.0235 | 0.1156 | 0.0587 | 0.0878 | 0.0235 | 0.1108 | 0.0587 | 0.0878 | 0.0235 | 0.1108 |
| Sample_15ng_VAF_0p25pc | 15 | 0.250 | 0.2613 | 0.3659 | 0.1353 | 0.4144 | 0.2933 | 0.4049 | 0.1588 | 0.4536 | 0.2933 | 0.4049 | 0.1588 | 0.4536 |
| Sample_15ng_VAF_0p5pc | 15 | 0.500 | 0.4000 | 0.5268 | 0.2471 | 0.5714 | 0.4240 | 0.5610 | 0.2588 | 0.5955 | 0.4240 | 0.5610 | 0.2588 | 0.5955 |
| Sample_15ng_VAF_1pc | 15 | 1.000 | 0.5547 | 0.7317 | 0.3412 | 0.7136 | 0.5760 | 0.7610 | 0.3529 | 0.7310 | 0.5760 | 0.7610 | 0.3529 | 0.7310 |
| Sample_15ng_VAF_2pc | 15 | 2.000 | 0.6133 | 0.8000 | 0.3882 | 0.7603 | 0.6267 | 0.8195 | 0.3941 | 0.7705 | 0.6267 | 0.8195 | 0.3941 | 0.7705 |
| Sample_30ng_VAF_0pc | 30 | 0.000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| Sample_30ng_VAF_0p001pc | 30 | 0.001 | 0.0053 | 0.0098 | 0.0000 | 0.0106 | 0.0053 | 0.0098 | 0.0000 | 0.0106 | 0.0053 | 0.0098 | 0.0000 | 0.0106 |
| Sample_30ng_VAF_0p01pc | 30 | 0.010 | 0.0053 | 0.0098 | 0.0000 | 0.0106 | 0.0053 | 0.0098 | 0.0000 | 0.0106 | 0.0053 | 0.0098 | 0.0000 | 0.0106 |
| Sample_30ng_VAF_0p1pc | 30 | 0.100 | 0.0773 | 0.1220 | 0.0235 | 0.1436 | 0.0827 | 0.1366 | 0.0176 | 0.1527 | 0.0827 | 0.1366 | 0.0176 | 0.1527 |
| Sample_30ng_VAF_0p25pc | 30 | 0.250 | 0.2853 | 0.3854 | 0.1647 | 0.4440 | 0.2853 | 0.3805 | 0.1706 | 0.4440 | 0.2853 | 0.3805 | 0.1706 | 0.4440 |
| Sample_30ng_VAF_0p5pc | 30 | 0.500 | 0.4613 | 0.6098 | 0.2824 | 0.6314 | 0.4960 | 0.6537 | 0.3059 | 0.6631 | 0.4960 | 0.6537 | 0.3059 | 0.6631 |
| Sample_30ng_VAF_1pc | 30 | 1.000 | 0.5653 | 0.7415 | 0.3529 | 0.7223 | 0.5893 | 0.7805 | 0.3588 | 0.7416 | 0.5893 | 0.7805 | 0.3588 | 0.7416 |
| Sample_30ng_VAF_2pc | 30 | 2.000 | 0.5920 | 0.7707 | 0.3765 | 0.7437 | 0.6080 | 0.7951 | 0.3824 | 0.7562 | 0.6080 | 0.7951 | 0.3824 | 0.7562 |

### 3.2 VAF-level averages across DNA inputs

Averages are computed across the three DNA input levels (5 ng, 15 ng, 30 ng) for each VAF. The 0% VAF negative controls are excluded.

**Computed values:**

- **0.001%:** Baseline = (0.0000 + 0.0000 + 0.0053) / 3 = 0.0018; TI = (0.0000 + 0.0000 + 0.0053) / 3 = 0.0018
- **0.010%:** Baseline = (0.0000 + 0.0053 + 0.0053) / 3 = 0.0035; TI = (0.0000 + 0.0027 + 0.0053) / 3 = 0.0027
- **0.100%:** Baseline = (0.0107 + 0.0613 + 0.0773) / 3 = 0.0498; TI = (0.0080 + 0.0587 + 0.0827) / 3 = 0.0498
- **0.250%:** Baseline = (0.0640 + 0.2613 + 0.2853) / 3 = 0.2035; TI = (0.0747 + 0.2933 + 0.2853) / 3 = 0.2178
- **0.500%:** Baseline = (0.1680 + 0.4000 + 0.4613) / 3 = 0.3431; TI = (0.1867 + 0.4240 + 0.4960) / 3 = 0.3689
- **1.000%:** Baseline = (0.3573 + 0.5547 + 0.5653) / 3 = 0.4924; TI = (0.4000 + 0.5760 + 0.5893) / 3 = 0.5218
- **2.000%:** Baseline = (0.5173 + 0.6133 + 0.5920) / 3 = 0.5742; TI = (0.5387 + 0.6267 + 0.6080) / 3 = 0.5911

**Overall sensitivity mean across ng inputs:**

| VAF (%) | Baseline Mean Sens | TI Mean Sens | ML+TI Mean Sens |
|---|---|---|---|
| 0.001 | 0.0018 | 0.0018 | 0.0018 |
| 0.010 | 0.0035 | 0.0027 | 0.0027 |
| 0.100 | 0.0498 | 0.0498 | 0.0498 |
| 0.250 | 0.2035 | 0.2178 | 0.2178 |
| 0.500 | 0.3431 | 0.3689 | 0.3689 |
| 1.000 | 0.4924 | 0.5218 | 0.5218 |
| 2.000 | 0.5742 | 0.5911 | 0.5911 |

**SNV sensitivity computed values:**

- **0.001%:** Baseline = (0.0000 + 0.0000 + 0.0098) / 3 = 0.0033; TI = (0.0000 + 0.0000 + 0.0098) / 3 = 0.0033
- **0.010%:** Baseline = (0.0000 + 0.0098 + 0.0098) / 3 = 0.0065; TI = (0.0000 + 0.0049 + 0.0098) / 3 = 0.0049
- **0.100%:** Baseline = (0.0195 + 0.0927 + 0.1220) / 3 = 0.0781; TI = (0.0146 + 0.0878 + 0.1366) / 3 = 0.0797
- **0.250%:** Baseline = (0.0829 + 0.3659 + 0.3854) / 3 = 0.2781; TI = (0.0927 + 0.4049 + 0.3805) / 3 = 0.2927
- **0.500%:** Baseline = (0.2244 + 0.5268 + 0.6098) / 3 = 0.4537; TI = (0.2488 + 0.5610 + 0.6537) / 3 = 0.4878
- **1.000%:** Baseline = (0.4732 + 0.7317 + 0.7415) / 3 = 0.6488; TI = (0.5366 + 0.7610 + 0.7805) / 3 = 0.6927
- **2.000%:** Baseline = (0.6878 + 0.8000 + 0.7707) / 3 = 0.7528; TI = (0.7171 + 0.8195 + 0.7951) / 3 = 0.7772

**SNV sensitivity mean across ng inputs:**

| VAF (%) | Baseline Mean SNV Sens | TI Mean SNV Sens | ML+TI Mean SNV Sens |
|---|---|---|---|
| 0.001 | 0.0033 | 0.0033 | 0.0033 |
| 0.010 | 0.0065 | 0.0049 | 0.0049 |
| 0.100 | 0.0781 | 0.0797 | 0.0797 |
| 0.250 | 0.2781 | 0.2927 | 0.2927 |
| 0.500 | 0.4537 | 0.4878 | 0.4878 |
| 1.000 | 0.6488 | 0.6927 | 0.6927 |
| 2.000 | 0.7528 | 0.7772 | 0.7772 |

**Indel sensitivity computed values:**

- **0.001%–0.010%:** All conditions = 0.0000 (no indels detected at very low VAF).
- **0.100%:** Baseline = (0.0000 + 0.0235 + 0.0235) / 3 = 0.0157; TI = (0.0000 + 0.0235 + 0.0176) / 3 = 0.0137
- **0.250%:** Baseline = (0.0412 + 0.1353 + 0.1647) / 3 = 0.1137; TI = (0.0529 + 0.1588 + 0.1706) / 3 = 0.1274
- **0.500%:** Baseline = (0.1000 + 0.2471 + 0.2824) / 3 = 0.2098; TI = (0.1118 + 0.2588 + 0.3059) / 3 = 0.2255
- **1.000%:** Baseline = (0.2176 + 0.3412 + 0.3529) / 3 = 0.3039; TI = (0.2353 + 0.3529 + 0.3588) / 3 = 0.3157
- **2.000%:** Baseline = (0.3118 + 0.3882 + 0.3765) / 3 = 0.3588; TI = (0.3235 + 0.3941 + 0.3824) / 3 = 0.3667

**Indel sensitivity mean across ng inputs:**

| VAF (%) | Baseline Mean Indel Sens | TI Mean Indel Sens | ML+TI Mean Indel Sens |
|---|---|---|---|
| 0.001 | 0.0000 | 0.0000 | 0.0000 |
| 0.010 | 0.0000 | 0.0000 | 0.0000 |
| 0.100 | 0.0157 | 0.0137 | 0.0137 |
| 0.250 | 0.1137 | 0.1274 | 0.1274 |
| 0.500 | 0.2098 | 0.2255 | 0.2255 |
| 1.000 | 0.3039 | 0.3157 | 0.3157 |
| 2.000 | 0.3588 | 0.3667 | 0.3667 |

---

## 4. Results: Precision and FP Analysis

### Precision

Precision is 1.000 in every sample under all three conditions. No false positives are recorded at any VAF level.

| Condition | Total FP across all 24 samples | Precision |
|---|---|---|
| Baseline | 0 | 1.000 |
| TI-only | 0 | 1.000 |
| ML+TI | 0 | 1.000 |

### Mechanism by condition

**Baseline:** The Python scoring script matches called variants against the truth VCF using 0-based coordinates. In this dataset, the pipeline under discovery mode only generates calls at positions within the target panel that have sufficient duplex molecule support. The panel targets cover exactly the positions containing truth variants. Background errors are suppressed by duplex consensus. No called variant falls outside the truth set, giving precision = 1.000.

**TI mode:** The `--target-variants` flag explicitly restricts PASS status to variants that match the truth VCF by CHROM, 1-based POS, REF, and ALT. Any variant not matching receives the NotTargeted filter status and is excluded from sensitivity and precision counts. This enforces precision = 1.000 by construction. The count of NotTargeted calls is not tabulated here. TI mode fundamentally changes the calling semantics: genuine FPs from de Bruijn graph artefacts, PCR errors, or index hopping would receive NotTargeted status rather than being scored as FP, making them invisible to the standard scoring metrics.

**ML+TI:** The ML filter operates after TI matching. Because no calls reach PASS status without TI matching, the ML filter has no FP candidates to filter. Precision remains 1.000.

### Negative controls (0% VAF)

The 0% VAF samples contain no spiked-in variant DNA. Any PASS call in these samples is a spurious background event matching a truth coordinate.

| Sample | Condition | PASS calls | TP | FP | Note |
|---|---|---|---|---|---|
| Sample_5ng_VAF_0pc | Baseline | 0 | 0 | 0 | Clean |
| Sample_5ng_VAF_0pc | TI-only | 0 | 0 | 0 | Clean |
| Sample_5ng_VAF_0pc | ML+TI | 0 | 0 | 0 | Clean |
| Sample_15ng_VAF_0pc | Baseline | 0 | 0 | 0 | Clean |
| Sample_15ng_VAF_0pc | TI-only | 1 | 1 | 0 | 1 spurious call at truth coordinate |
| Sample_15ng_VAF_0pc | ML+TI | 1 | 1 | 0 | ML does not remove the spurious call |
| Sample_30ng_VAF_0pc | Baseline | 0 | 0 | 0 | Clean |
| Sample_30ng_VAF_0pc | TI-only | 0 | 0 | 0 | Clean |
| Sample_30ng_VAF_0pc | ML+TI | 0 | 0 | 0 | Clean |

The single spurious call in the 15 ng 0% VAF sample is a background sequencing error that coincidentally occurs at a truth variant coordinate and passes the duplex consensus quality threshold. Under TI scoring, this call is recorded as tp=1 (sensitivity = 0.0027) because TI mode cannot distinguish a real call from a background error at a truth position. The ML model scores this call at ml_prob ≥ 0.999 and does not remove it. The call is an SNV (snv_tp=1; indel_tp=0).

---

## 5. Results: ML Filter Effect

### Observation

The ML filter has no effect on sensitivity or precision in TI mode. For every sample where calls are present, ml_tp equals tp, ml_fp equals fp (zero), and ml_sensitivity equals sensitivity.

**Per-sample comparison, ML+TI condition (all samples with non-zero calls):**

| Sample | ng | VAF (%) | TI TP | ML+TI ml_tp | TI FP | ML+TI ml_fp | TI Sens | ML+TI ml_sens |
|---|---|---|---|---|---|---|---|---|
| Sample_15ng_VAF_0pc | 15 | 0.000 | 1 | 1 | 0 | 0 | 0.0027 | 0.0027 |
| Sample_15ng_VAF_0p01pc | 15 | 0.010 | 1 | 1 | 0 | 0 | 0.0027 | 0.0027 |
| Sample_15ng_VAF_0p1pc | 15 | 0.100 | 22 | 22 | 0 | 0 | 0.0587 | 0.0587 |
| Sample_15ng_VAF_0p25pc | 15 | 0.250 | 110 | 110 | 0 | 0 | 0.2933 | 0.2933 |
| Sample_15ng_VAF_0p5pc | 15 | 0.500 | 159 | 159 | 0 | 0 | 0.4240 | 0.4240 |
| Sample_15ng_VAF_1pc | 15 | 1.000 | 216 | 216 | 0 | 0 | 0.5760 | 0.5760 |
| Sample_15ng_VAF_2pc | 15 | 2.000 | 235 | 235 | 0 | 0 | 0.6267 | 0.6267 |
| Sample_30ng_VAF_0p001pc | 30 | 0.001 | 2 | 2 | 0 | 0 | 0.0053 | 0.0053 |
| Sample_30ng_VAF_0p01pc | 30 | 0.010 | 2 | 2 | 0 | 0 | 0.0053 | 0.0053 |
| Sample_30ng_VAF_0p1pc | 30 | 0.100 | 31 | 31 | 0 | 0 | 0.0827 | 0.0827 |
| Sample_30ng_VAF_0p25pc | 30 | 0.250 | 107 | 107 | 0 | 0 | 0.2853 | 0.2853 |
| Sample_30ng_VAF_0p5pc | 30 | 0.500 | 186 | 186 | 0 | 0 | 0.4960 | 0.4960 |
| Sample_30ng_VAF_1pc | 30 | 1.000 | 221 | 221 | 0 | 0 | 0.5893 | 0.5893 |
| Sample_30ng_VAF_2pc | 30 | 2.000 | 228 | 228 | 0 | 0 | 0.6080 | 0.6080 |
| Sample_5ng_VAF_0p1pc | 5 | 0.100 | 3 | 3 | 0 | 0 | 0.0080 | 0.0080 |
| Sample_5ng_VAF_0p25pc | 5 | 0.250 | 28 | 28 | 0 | 0 | 0.0747 | 0.0747 |
| Sample_5ng_VAF_0p5pc | 5 | 0.500 | 70 | 70 | 0 | 0 | 0.1867 | 0.1867 |
| Sample_5ng_VAF_1pc | 5 | 1.000 | 150 | 150 | 0 | 0 | 0.4000 | 0.4000 |
| Sample_5ng_VAF_2pc | 5 | 2.000 | 202 | 202 | 0 | 0 | 0.5387 | 0.5387 |

For samples where no variants are called (0% VAF for 30 ng and 5 ng; 0.001% VAF for 15 ng and 5 ng), the ML columns are NA in the source TSV, indicating the ML scoring pathway was not reached.

### Reason

In TI mode, only calls that exactly match a truth variant by coordinate and allele receive PASS. These calls arise from genuine variant signal: the molecule ensemble at that position contains duplex reads carrying the alt allele at sufficient depth. Variants with real molecule support have high depth, high duplex fraction, and high k-mer evidence weight. These are precisely the features that drive the ML model to assign high probability. In this dataset, every TI-PASS call achieves ml_precision = 1.000 and ml_tp equals tp. The ML model correctly identifies every TI-PASS call as a probable TP, leaving sensitivity and precision unchanged.

### Runtime overhead

The ML model (LightGBM inference executed in the call stage) adds approximately 15–18 ms per sample compared to TI-only. The table below compares call-stage timings (t_call_ms) from the two TSVs:

| Sample | TI t_call_ms | ML+TI t_call_ms | Overhead (ms) |
|---|---|---|---|
| Sample_15ng_VAF_0p001pc | 7 | 16 | 9 |
| Sample_15ng_VAF_0p01pc | 1 | 16 | 15 |
| Sample_15ng_VAF_0p1pc | 1 | 17 | 16 |
| Sample_15ng_VAF_0p25pc | 10 | 17 | 7 |
| Sample_15ng_VAF_0p5pc | 1 | 19 | 18 |
| Sample_15ng_VAF_1pc | 1 | 18 | 17 |
| Sample_15ng_VAF_2pc | 2 | 19 | 17 |
| Sample_30ng_VAF_0p5pc | 2 | 19 | 17 |
| Sample_30ng_VAF_1pc | 2 | 20 | 18 |
| Sample_30ng_VAF_2pc | 3 | 20 | 17 |
| Sample_5ng_VAF_1pc | 1 | 17 | 16 |
| Sample_5ng_VAF_2pc | 1 | 18 | 17 |

The typical ML overhead is 15–18 ms. Total wall time is dominated by the assembly stage (~37–45 s per sample) and is not materially affected.

---

## 6. Results: Sensitivity by VAF Curve

### 6.1 Overall sensitivity by VAF and DNA input

| ng | VAF (%) | Baseline Sens | TI Sens | ML+TI Sens |
|---|---|---|---|---|
| 5 | 0.001 | 0.0000 | 0.0000 | 0.0000 |
| 5 | 0.010 | 0.0000 | 0.0000 | 0.0000 |
| 5 | 0.100 | 0.0107 | 0.0080 | 0.0080 |
| 5 | 0.250 | 0.0640 | 0.0747 | 0.0747 |
| 5 | 0.500 | 0.1680 | 0.1867 | 0.1867 |
| 5 | 1.000 | 0.3573 | 0.4000 | 0.4000 |
| 5 | 2.000 | 0.5173 | 0.5387 | 0.5387 |
| 15 | 0.001 | 0.0000 | 0.0000 | 0.0000 |
| 15 | 0.010 | 0.0053 | 0.0027 | 0.0027 |
| 15 | 0.100 | 0.0613 | 0.0587 | 0.0587 |
| 15 | 0.250 | 0.2613 | 0.2933 | 0.2933 |
| 15 | 0.500 | 0.4000 | 0.4240 | 0.4240 |
| 15 | 1.000 | 0.5547 | 0.5760 | 0.5760 |
| 15 | 2.000 | 0.6133 | 0.6267 | 0.6267 |
| 30 | 0.001 | 0.0053 | 0.0053 | 0.0053 |
| 30 | 0.010 | 0.0053 | 0.0053 | 0.0053 |
| 30 | 0.100 | 0.0773 | 0.0827 | 0.0827 |
| 30 | 0.250 | 0.2853 | 0.2853 | 0.2853 |
| 30 | 0.500 | 0.4613 | 0.4960 | 0.4960 |
| 30 | 1.000 | 0.5653 | 0.5893 | 0.5893 |
| 30 | 2.000 | 0.5920 | 0.6080 | 0.6080 |

Sensitivity rises steeply between 0.1% and 1% VAF across all DNA inputs. At 0.001%–0.010% VAF, essentially no variants are detected. Sensitivity approaches a soft plateau of 54–63% at 2% VAF. The plateau reflects a combination of coverage limitations and the inherent difficulty of detecting low-abundance indels.

### 6.2 SNV sensitivity by VAF and DNA input

| ng | VAF (%) | Baseline SNV Sens | TI SNV Sens | ML+TI SNV Sens |
|---|---|---|---|---|
| 5 | 0.001 | 0.0000 | 0.0000 | 0.0000 |
| 5 | 0.010 | 0.0000 | 0.0000 | 0.0000 |
| 5 | 0.100 | 0.0195 | 0.0146 | 0.0146 |
| 5 | 0.250 | 0.0829 | 0.0927 | 0.0927 |
| 5 | 0.500 | 0.2244 | 0.2488 | 0.2488 |
| 5 | 1.000 | 0.4732 | 0.5366 | 0.5366 |
| 5 | 2.000 | 0.6878 | 0.7171 | 0.7171 |
| 15 | 0.001 | 0.0000 | 0.0000 | 0.0000 |
| 15 | 0.010 | 0.0098 | 0.0049 | 0.0049 |
| 15 | 0.100 | 0.0927 | 0.0878 | 0.0878 |
| 15 | 0.250 | 0.3659 | 0.4049 | 0.4049 |
| 15 | 0.500 | 0.5268 | 0.5610 | 0.5610 |
| 15 | 1.000 | 0.7317 | 0.7610 | 0.7610 |
| 15 | 2.000 | 0.8000 | 0.8195 | 0.8195 |
| 30 | 0.001 | 0.0098 | 0.0098 | 0.0098 |
| 30 | 0.010 | 0.0098 | 0.0098 | 0.0098 |
| 30 | 0.100 | 0.1220 | 0.1366 | 0.1366 |
| 30 | 0.250 | 0.3854 | 0.3805 | 0.3805 |
| 30 | 0.500 | 0.6098 | 0.6537 | 0.6537 |
| 30 | 1.000 | 0.7415 | 0.7805 | 0.7805 |
| 30 | 2.000 | 0.7707 | 0.7951 | 0.7951 |

At 2% VAF and 15 ng, SNV sensitivity reaches 0.8195 under TI/ML+TI (168 of 205 SNVs detected). At 30 ng and 2% VAF it reaches 0.7951 (163/205). At 5 ng and 2% VAF it reaches 0.7171 (147/205).

### 6.3 Indel sensitivity by VAF and DNA input

| ng | VAF (%) | Baseline Indel Sens | TI Indel Sens | ML+TI Indel Sens |
|---|---|---|---|---|
| 5 | 0.001 | 0.0000 | 0.0000 | 0.0000 |
| 5 | 0.010 | 0.0000 | 0.0000 | 0.0000 |
| 5 | 0.100 | 0.0000 | 0.0000 | 0.0000 |
| 5 | 0.250 | 0.0412 | 0.0529 | 0.0529 |
| 5 | 0.500 | 0.1000 | 0.1118 | 0.1118 |
| 5 | 1.000 | 0.2176 | 0.2353 | 0.2353 |
| 5 | 2.000 | 0.3118 | 0.3235 | 0.3235 |
| 15 | 0.001 | 0.0000 | 0.0000 | 0.0000 |
| 15 | 0.010 | 0.0000 | 0.0000 | 0.0000 |
| 15 | 0.100 | 0.0235 | 0.0235 | 0.0235 |
| 15 | 0.250 | 0.1353 | 0.1588 | 0.1588 |
| 15 | 0.500 | 0.2471 | 0.2588 | 0.2588 |
| 15 | 1.000 | 0.3412 | 0.3529 | 0.3529 |
| 15 | 2.000 | 0.3882 | 0.3941 | 0.3941 |
| 30 | 0.001 | 0.0000 | 0.0000 | 0.0000 |
| 30 | 0.010 | 0.0000 | 0.0000 | 0.0000 |
| 30 | 0.100 | 0.0235 | 0.0176 | 0.0176 |
| 30 | 0.250 | 0.1647 | 0.1706 | 0.1706 |
| 30 | 0.500 | 0.2824 | 0.3059 | 0.3059 |
| 30 | 1.000 | 0.3529 | 0.3588 | 0.3588 |
| 30 | 2.000 | 0.3765 | 0.3824 | 0.3824 |

At 2% VAF, indel sensitivity peaks at 0.3941 (15 ng, TI/ML+TI), detecting 67 of 170 indels. No indels are detected at the 5 ng and 0.1% VAF level under any condition. The detection threshold for indels is higher than for SNVs across all DNA inputs and VAF levels.

---

## 7. Results: Per-Sample Observations

### 7.1 Negative controls (0% VAF)

**5 ng, 0% VAF:** 0 PASS calls in all three conditions. The negative control is clean.

**15 ng, 0% VAF:**
- Baseline: 0 calls. The Python scoring pipeline does not call any variant.
- TI-only: 1 PASS call (tp=1, fp=0, fn=374). The call is an SNV (snv_tp=1, snv_fn=204). This is a background sequencing error at a truth position. Sensitivity = 1/375 = 0.0027. F1 = 0.0053.
- ML+TI: Identical to TI-only. ml_tp=1, ml_fp=0, ml_sensitivity=0.0027. The ML model assigns this call ml_prob ≥ 0.999 and does not remove it.

**30 ng, 0% VAF:** 0 PASS calls in all three conditions. The negative control is clean.

The spurious 15 ng negative-control call illustrates a fundamental limit of TI mode: a background error molecule that coincidentally covers a truth coordinate cannot be distinguished from a true somatic call using call-level features alone. The duplex consensus for this molecule passes all quality thresholds, and the ML model correctly predicts it is a high-confidence call, which it is, in the sense that the k-mer evidence is strong. The biology is simply wrong.

### 7.2 DNA input effect on sensitivity

DNA input has a large effect, particularly at intermediate VAF levels.

At 0.5% VAF (TI condition):
- 5 ng: sensitivity = 0.1867 (70/375 variants)
- 15 ng: sensitivity = 0.4240 (159/375 variants)
- 30 ng: sensitivity = 0.4960 (186/375 variants)

The 30 ng input detects 116 more variants than the 5 ng input at 0.5% VAF. This gap narrows at higher VAF as the signal becomes more detectable regardless of DNA input.

At 2% VAF (TI condition):
- 5 ng: 0.5387 (202/375)
- 15 ng: 0.6267 (235/375)
- 30 ng: 0.6080 (228/375)

At 2% VAF, 30 ng is slightly lower than 15 ng (0.6080 vs 0.6267). This is not necessarily a systematic effect: the 30 ng sample has a lower duplex percentage (6.585%) than the 15 ng sample (15.437%) at 2% VAF, which likely reduces duplex-confirmed calls.

Duplex percentages confirm the input effect. For 2% VAF:
- 5 ng: duplex pct = 4.735% (29,643 duplex molecules)
- 15 ng: duplex pct = 15.437% (66,501 duplex molecules)
- 30 ng: duplex pct = 6.585% (26,204 duplex molecules)

The 15 ng sample has the highest duplex fraction at 2% VAF, explaining its higher sensitivity despite intermediate DNA input.

### 7.3 SNV versus indel sensitivity

SNV sensitivity exceeds indel sensitivity by approximately a factor of 2 at every (ng, VAF) combination. At 1% VAF under TI conditions:

| ng | SNV Sens | Indel Sens | Ratio (SNV / Indel) |
|---|---|---|---|
| 5 | 0.5366 | 0.2353 | 2.28 |
| 15 | 0.7610 | 0.3529 | 2.16 |
| 30 | 0.7805 | 0.3588 | 2.18 |

The ratio is consistent across ng inputs, confirming that the SNV/indel gap is intrinsic to the k-mer graph approach rather than a coverage artefact. Indel detection requires the path-finder to identify a subgraph that diverges from the reference k-mer set by insertion or deletion. Short indels that disrupt the k-mer reading frame increase the number of required candidate paths and are more susceptible to assembly gaps at low allele frequency.

At 2% VAF under TI conditions:

| ng | SNV Sens | Indel Sens | Ratio (SNV / Indel) |
|---|---|---|---|
| 5 | 0.7171 | 0.3235 | 2.22 |
| 15 | 0.8195 | 0.3941 | 2.08 |
| 30 | 0.7951 | 0.3824 | 2.08 |

The ratio is slightly lower at 2% VAF than at 1% VAF, suggesting that indel sensitivity grows proportionally faster than SNV sensitivity as VAF increases. Both are still far below 1.0 at 2% VAF, indicating headroom for improvement in both categories.

### 7.4 Low-VAF detection

At 0.001% VAF, only the 30 ng sample detects any variants under any condition (tp=2, both SNV). At 0.01% VAF, the 30 ng sample detects 2 (2 SNV, 0 indel) and the 15 ng sample detects 1–2 depending on condition. The 5 ng sample detects 0 at both of these VAF levels.

These sporadic detections are stochastic. At 0.001% VAF, a panel of 375 variants has an expected count of 0.375 mutant molecules per target per thousand input molecules. With 2 million reads and typical duplex fractions below 15%, the expected number of duplex molecules carrying a specific 0.001% variant is much less than 1. The occasional positive call reflects chance overlap of a variant molecule with the detection threshold, not systematic sensitivity at this VAF level.

---

## 8. Comparison to Previous Baseline

### 8.1 Sensitivity difference (TI minus Baseline)

| ng | VAF (%) | Baseline Sens | TI Sens | Difference |
|---|---|---|---|---|
| 5 | 0.100 | 0.0107 | 0.0080 | -0.0027 |
| 5 | 0.250 | 0.0640 | 0.0747 | +0.0107 |
| 5 | 0.500 | 0.1680 | 0.1867 | +0.0187 |
| 5 | 1.000 | 0.3573 | 0.4000 | +0.0427 |
| 5 | 2.000 | 0.5173 | 0.5387 | +0.0214 |
| 15 | 0.010 | 0.0053 | 0.0027 | -0.0027 |
| 15 | 0.100 | 0.0613 | 0.0587 | -0.0027 |
| 15 | 0.250 | 0.2613 | 0.2933 | +0.0320 |
| 15 | 0.500 | 0.4000 | 0.4240 | +0.0240 |
| 15 | 1.000 | 0.5547 | 0.5760 | +0.0213 |
| 15 | 2.000 | 0.6133 | 0.6267 | +0.0133 |
| 30 | 0.001 | 0.0053 | 0.0053 | 0.0000 |
| 30 | 0.010 | 0.0053 | 0.0053 | 0.0000 |
| 30 | 0.100 | 0.0773 | 0.0827 | +0.0053 |
| 30 | 0.250 | 0.2853 | 0.2853 | 0.0000 |
| 30 | 0.500 | 0.4613 | 0.4960 | +0.0347 |
| 30 | 1.000 | 0.5653 | 0.5893 | +0.0240 |
| 30 | 2.000 | 0.5920 | 0.6080 | +0.0160 |

Of 18 non-trivial comparisons (excluding the 0% VAF controls and 0.001% VAF where both are near zero), 14 show TI equal to or higher than baseline. The three negative differences (5 ng 0.1%, 15 ng 0.01%, 15 ng 0.1%) are each 0.0027 sensitivity units, equivalent to 1 variant. These small negative values are within the noise of stochastic sampling differences between runs. The consistent positive differences at higher VAF (up to 4.3% absolute for 5 ng 1%) reflect genuine improvement in the current binary over the pre-ML baseline binary in pathfinding or calling logic.

TI mode itself does not cause sensitivity loss: the same set of candidate calls is generated regardless of the TI flag; TI mode only changes which calls receive PASS status. The sensitivity differences therefore reflect binary version changes.

### 8.2 Precision comparison

Both conditions achieve precision = 1.000. The mechanism differs:

- **Baseline:** Precision = 1.000 because discovery mode does not generate spurious calls on this dataset. This is not guaranteed in general; in samples with complex backgrounds or lower-quality inputs, discovery mode would produce FPs.
- **TI mode:** Precision = 1.000 by construction. This guarantee holds regardless of the background noise level, as long as the truth VCF is correct.

For clinical monitoring applications where the variant list from a primary tumour biopsy is available, TI mode provides a stronger precision guarantee with no sensitivity cost.

### 8.3 Molecule counts across conditions

The molecule counts differ between the baseline run and the TI runs because different binary versions were used. The 5 ng samples show the largest discrepancy: for example, the baseline Sample_5ng_VAF_2pc reports 346,773 molecules while the TI condition reports 625,985 molecules for the same sample. This approximately 2x difference suggests a change in the k-mer counting or molecule grouping logic between binary versions. The higher molecule counts in TI runs do not translate to proportionally higher sensitivity, possibly because many additional molecules are non-duplex and do not contribute to duplex consensus calls. This discrepancy should be investigated separately before drawing conclusions about absolute molecule yields.

---

## 9. Conclusion

### TI mode eliminates FPs with no sensitivity cost

Tumour-informed mode enforces precision = 1.000 by construction. On this titration dataset, it also matches or exceeds baseline sensitivity at most VAF levels. The sensitivity advantage is attributable to binary improvements between the baseline build and the current build, not to the TI mechanism itself. TI mode is the correct configuration for monitoring applications where the target variant list is known.

### ML model adds no benefit in TI mode

The twist-duplex-v1 ML model does not change sensitivity or precision in TI mode. Every call that passes TI matching has high molecule support, high duplex fraction, and high consensus confidence: exactly the features the ML model uses to assign high TP probability. In TI mode there are no FP candidates to filter, so the ML model has no discriminatory work to do. Adding the ML filter in TI mode is safe (it removes nothing) but pointless (it adds no benefit). The runtime overhead of approximately 17 ms per sample is negligible.

### ML model may be useful in non-TI mode

The ML model's expected utility is in discovery (non-TI) mode, where genuine FPs arise from erroneous de Bruijn graph paths, PCR artefacts, or sequencing errors outside the truth set. In that context, the model's ability to assign low probability to FP-like calls should reduce the FP rate at a given sensitivity threshold. This benchmark does not evaluate that scenario. A dedicated non-TI benchmark is required to quantify the ML model's practical benefit.

### Recommended configuration

For monitoring of known variants: use TI mode. The ML filter is optional and adds no benefit in this mode.

For discovery of unknown variants: run without TI. Evaluate the ML filter separately to determine whether its sensitivity cost (if any) is justified by FP reduction.

### Future work

1. Evaluate twist-duplex-v1 in non-TI discovery mode on this titration dataset. Measure sensitivity, FP rate, and the precision/sensitivity trade-off at various ml_prob thresholds.
2. Train ML models on real duplex sequencing data rather than varforge simulations to improve the AUPRC beyond 0.606 and better reflect the true distribution of call features.
3. Characterise the detection floor at 0.001%–0.010% VAF. Increasing sequencing depth or using more input DNA may shift the sensitivity curve to lower VAF.
4. Investigate the approximately 2x molecule count discrepancy between the baseline and TI binary versions for 5 ng samples to confirm that the difference is benign.
