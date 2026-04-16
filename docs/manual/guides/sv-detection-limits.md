# SV Detection Limits

## Overview

Detection limits define the lowest VAF at which kam reliably identifies each SV type. Knowing
these limits is critical for clinical applications: a monitoring assay must know whether a
negative result means the variant is absent or simply below the detection floor.

This guide covers how to measure detection limits using synthetic data and how to interpret the
results.

---

## What detection limits mean clinically

A detection limit is the minimum VAF at which a given SV type is detected with at least 80%
sensitivity. Below this threshold, the assay may miss true positives at an unacceptable rate.

For ctDNA monitoring, the detection limit determines:

- Whether a rising tumour fraction is caught early enough for clinical intervention.
- Whether molecular residual disease (MRD) can be tracked at post-treatment levels.
- The confidence you can place in a negative result.

---

## Current approximate limits

These values come from synthetic benchmarks at 5000x coverage with Twist duplex UMI chemistry.

| SV type | Discovery mode | Tumour-informed mode |
|---------|---------------|---------------------|
| LargeDeletion | 0.10% VAF | 0.05% VAF |
| TandemDuplication | 0.25% VAF | 0.10% VAF |
| Inversion | 0.10% VAF | 0.05% VAF |
| NovelInsertion | 0.25% VAF | 0.10% VAF |
| InvDel | 0.25% VAF | 0.15% VAF |
| Fusion | 0.10% VAF | 0.05% VAF |

Tumour-informed mode achieves lower detection limits because it does not require the variant
to pass confidence filters independently. The caller already knows the expected allele and
applies a lower evidence threshold.

Detection limits improve with higher sequencing depth. At 10,000x coverage, limits drop by
roughly 2-fold.

---

## Running ultra-low VAF benchmarks

The standard SV benchmark covers 0.05% to 10% VAF. To measure detection limits more precisely,
run the ultra-low VAF extension which covers 0.01% to 0.04%.

### Step 1: generate ultra-low VAF data

```bash
python3 docs/benchmarking/03-sv-extended/scripts/make_invdel_suite.py \
  --vaf-levels 0.01,0.02,0.03,0.04 \
  --coverage 5000 \
  --output-dir bigdata/benchmarking/sv_ultravaf/configs/
```

Repeat for each SV type suite script (`make_novins_suite.py`, `make_fusion_suite.py`).

### Step 2: simulate reads

```bash
bash docs/benchmarking/03-sv-extended/scripts/run_sv_new_suite.sh \
  --config-dir bigdata/benchmarking/sv_ultravaf/configs/ \
  --output-dir bigdata/benchmarking/sv_ultravaf/results/
```

### Step 3: score results

```bash
python3 docs/benchmarking/03-sv-extended/scripts/score_sv_new_suite.py \
  --results-dir bigdata/benchmarking/sv_ultravaf/results/ \
  --output bigdata/benchmarking/sv_ultravaf/detection_limits.tsv
```

---

## Reading detection limit curves

The scoring script produces a sensitivity-vs-VAF curve for each SV type. The x-axis is VAF
(log scale, 0.01% to 10%). The y-axis is sensitivity (0% to 100%).

### Key features to look for

**The 80% sensitivity crossing.** The VAF at which the curve crosses 80% is the detection
limit. Values below this are unreliable for clinical use.

**The 100% sensitivity plateau.** The VAF above which every true positive is detected. For
clinical reporting, this is the "safe zone" where a negative result means the variant is
genuinely absent.

**Separation between discovery and tumour-informed.** The gap between the two curves shows the
benefit of prior knowledge. A large gap means tumour-informed mode significantly extends the
assay's reach.

**Type-specific differences.** Large deletions and inversions reach lower detection limits
than tandem duplications and novel insertions. This reflects the number of variant-specific
k-mers available for each type: deletions remove many k-mers, creating a strong signal.
Insertions add new k-mers that may overlap background sequences.

---

## Factors that affect detection limits

| Factor | Effect on limit |
|--------|----------------|
| Sequencing depth | Higher depth lowers the limit. 10,000x halves the limit compared to 5,000x. |
| Duplex fraction | Higher duplex fraction provides stronger per-molecule evidence. |
| Target window length | Longer targets provide more k-mers, improving sensitivity. |
| K-mer size | Smaller k increases sensitivity at ultra-low VAF but may increase false positives. |
| SV type | Deletions and inversions have lower limits than insertions and duplications. |
| Background error rate | Higher error rates raise the effective limit. Use `--max-vaf 0.35` to suppress germline noise. |

---

## Recommendations

- **Report detection limits per SV type**, not a single number for all SVs. The limits vary
  by type.
- **Use tumour-informed mode** for monitoring known SVs. It provides a 2-fold improvement in
  detection limits.
- **Run a 0% VAF control** alongside your detection limit experiment to measure the false
  positive baseline.
- **Consider the clinical question.** If the assay must detect MRD at 0.01% VAF, current SV
  detection cannot reliably do this. SNV detection with duplex confirmation is more sensitive
  at extreme dilutions.
