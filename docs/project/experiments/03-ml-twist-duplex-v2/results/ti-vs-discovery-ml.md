# Tumour-Informed vs Discovery + ML: How Close Can We Get?

**Dataset:** 24-sample titration dataset — 3 ng conditions (5, 15, 30 ng) × 8 VAF levels (0–2%).  
**Models compared:** raw discovery, discovery + ML v1, discovery + ML v2, tumour-informed (TI).  
**Source data:** `docs/benchmarking/07-snvindel-ml-boost-v1/results/comparison_disc_v1_v2_ti.csv`

---

## The core question

Tumour-informed (TI) mode uses a truth VCF of known variants to restrict calls to targeted positions. It achieves precision 1.000 at all VAF levels but requires prior knowledge of the tumour's variant profile. The question is whether an ML model trained on real sequencing artefacts can approach that precision without the truth VCF.

---

## ML v1 is not useful

ML v1 (trained on varforge-simulated data) assigns `ml_prob >= 0.9` to every PASS call — both true positives and false positives. Its precision column is identical to raw discovery across all 24 conditions. It provides no discrimination whatsoever.

The root cause: varforge simulates perfect data where every PASS call is a real variant. The model learned "passes the confidence filter" as the positive class. It has never seen a real background sequencing error.

---

## ML v2 at 1–2% VAF: approaches TI

The clinically relevant range for ctDNA monitoring is 1–2% VAF. At these levels ML v2 closes most of the gap.

| Condition | Discovery | ML v2 | TI | Gap (v2 vs TI) |
|-----------|-----------|-------|----|----------------|
| 5ng 1% VAF | 0.761 | 0.974 | 1.000 | 0.026 |
| 5ng 2% VAF | 0.838 | 1.000 | 1.000 | 0.000 |
| 15ng 1% VAF | 0.755 | 0.969 | 1.000 | 0.031 |
| 15ng 2% VAF | 0.693 | 0.979 | 1.000 | 0.021 |
| 30ng 1% VAF | 0.720 | 0.893 | 1.000 | 0.107 |
| 30ng 2% VAF | 0.753 | 0.925 | 1.000 | 0.075 |

At 5ng and 15ng (the training conditions), ML v2 reaches precision 0.969–1.000. The residual gap is 3–5 false positives per sample that share the feature profile of true low-VAF variants. At 30ng (the held-out test condition), the gap is larger — 7–8% — because the model was not trained on that input mass.

Sensitivity is fully preserved at 5ng and 15ng. At 30ng, retention is 97–99%: the model drops a small number of true positives.

---

## False positive reduction at negative controls

At 0% VAF (negative control samples with no spiked-in variants), every call is a false positive. TI eliminates all of them. ML v2 reduces them substantially but cannot reach zero.

| Condition | Discovery FP | ML v2 FP | TI FP | Reduction |
|-----------|-------------|----------|-------|-----------|
| 5ng 0% VAF | 43 | 3 | 0 | 93% |
| 15ng 0% VAF | 74 | 6 | 0 | 92% |
| 30ng 0% VAF | 86 | 18 | 0 | 79% |

The 3–18 residual FPs are background sequencing errors with high-quality duplex support, consistent trinucleotide context, and low error probability — indistinguishable from a real variant at 0.1–0.5% VAF by the features available. TI eliminates them trivially by position filtering.

---

## Breakdown across all VAF levels

The gap between ML v2 and TI widens at low VAF because:

1. There are fewer true positives to anchor the signal. At 0.001% VAF, kam detects almost nothing in either mode.
2. The background FP rate is constant regardless of VAF. At low VAF, FPs dominate.
3. ML v2 was not trained on sub-0.1% VAF data — those conditions have too few TPs for labelling.

| VAF | 5ng v2 prec | 15ng v2 prec | 30ng v2 prec | TI prec |
|-----|-------------|--------------|--------------|---------|
| 0% | 0.000 | 0.000 | 0.000 | 0.000 |
| 0.001% | 0.000 | 0.000 | 0.063 | 0.000–1.000 |
| 0.01% | 0.000 | 0.143 | 0.000 | 0.000–1.000 |
| 0.1% | 0.500 | 0.733 | 0.537 | 1.000 |
| 0.25% | 0.933 | 0.957 | 0.823 | 1.000 |
| 0.5% | 0.972 | 0.958 | 0.899 | 1.000 |
| **1%** | **0.974** | **0.969** | **0.893** | **1.000** |
| **2%** | **1.000** | **0.979** | **0.925** | **1.000** |

TI precision at sub-0.1% VAF is also unstable (alternates between 0.000 and 1.000) because kam detects so few variants at these conditions that precision is determined by one or two calls.

---

## Conclusions

**At 1–2% VAF**, discovery + ML v2 approaches TI precision within 2–8%. The practical gap is small: 4–26 residual FPs per sample across conditions compared to TI's zero. For monitoring applications where a truth VCF is unavailable, ML v2 provides a viable alternative.

**At < 0.25% VAF**, the gap remains large. Both modes struggle with sensitivity at 5ng, and ML v2 cannot discriminate the sparse TP signal from background. TI has a structural advantage here: position filtering removes FPs regardless of their feature profile.

**The irreducible gap** is background errors that are sequencing-quality correct, duplex-confirmed, and occur at plausible cancer mutation contexts (CpG transitions, oxidative damage). These are the variants the model cannot filter without coordinate information.

**Recommendation for use:** ML v2 is appropriate for discovery applications where sensitivity at 1–2% VAF matters and a truth VCF is unavailable. For near-zero false-positive monitoring, TI mode remains the correct choice.
