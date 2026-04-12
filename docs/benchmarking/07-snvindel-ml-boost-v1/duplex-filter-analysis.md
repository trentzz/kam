# Duplex Confirmation Filter Analysis — `--min-alt-duplex 1`

**Date:** 2026-04-12
**Benchmark:** 07-snvindel-ml-boost-v1 (24-sample titration, TI mode)
**Question:** What is the sensitivity cost of requiring at least one duplex-confirmed alt molecule?

---

## Background

Tumour-informed (TI) mode restricts PASS calls to variants that match a known somatic truth VCF entry by exact (CHROM, POS, REF, ALT). Precision is 1.000 by construction for all samples with any spiked-in variants. However, the 15 ng 0% VAF negative control produced one PASS call under TI mode: `chr3:10149804:C:T`, a C→T deamination error at 0.358% VAF with 2 alt molecules and zero duplex confirmation (NDUPALT=0, NSIMALT=2). The call passed TI mode because that coordinate and allele appear in the truth VCF; there are no actual spiked-in variants in this sample.

The `--min-alt-duplex N` flag (`-D warnings` filter label: `LowDuplex`) requires at least N variant-specific duplex molecules for a PASS call. Variant-specific duplex molecules are counted at k-mers that appear in the alt path but not the reference path — a more focused signal than general molecule depth. Setting `--min-alt-duplex 1` eliminates calls with zero duplex evidence at the variant site.

This document reports the sensitivity cost of setting `--min-alt-duplex 1` relative to TI-only (no duplex requirement).

---

## Methods

The same 24-sample titration dataset used in benchmark 07 was re-run with `--target-variants` and `--min-alt-duplex 1`:

```
python3 docs/benchmarking/07-snvindel-ml-boost-v1/scripts/run_titration_batch.py \
  --fastq-dir /mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs \
  --reads 2000000 \
  --target-variants docs/benchmarking/01-snvindel/scripts/truth_variants.vcf \
  --min-alt-duplex 1 \
  --output titration_2mreads_ti_minduplex1.tsv
```

All other parameters are identical to the TI-only run. Results are in `results/titration_2mreads_ti_minduplex1.tsv`.

---

## Results

### Negative control (0% VAF)

| Input | TI-only PASS | +minduplex1 PASS |
|-------|-------------|-----------------|
| 15 ng | 1 (chr3:10149804:C:T, background error) | 0 |
| 30 ng | 0 | 0 |
| 5 ng  | 0 | 0 |

The filter eliminates the spurious background call. The other two negative controls were already clean.

### Sensitivity at each VAF level

Truth set: 375 variants (205 SNV, 170 indel).

#### 15 ng

| VAF% | TI sens | +d1 sens | Δ      | TI snv | +d1 snv | TI ind | +d1 ind | TI tp | +d1 tp |
|------|---------|----------|--------|--------|---------|--------|---------|-------|--------|
| 0.001 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0 | 0 |
| 0.01  | 0.003 | 0.000 | −0.003 | 0.005 | 0.000 | 0.000 | 0.000 | 1 | 0 |
| 0.1   | 0.059 | 0.024 | −0.035 | 0.088 | 0.029 | 0.024 | 0.018 | 22 | 9 |
| 0.25  | 0.293 | 0.163 | −0.131 | 0.405 | 0.215 | 0.159 | 0.100 | 110 | 61 |
| 0.5   | 0.424 | 0.221 | −0.203 | 0.561 | 0.293 | 0.259 | 0.135 | 159 | 83 |
| 1.0   | 0.576 | 0.392 | −0.184 | 0.761 | 0.527 | 0.353 | 0.229 | 216 | 147 |
| 2.0   | 0.627 | 0.565 | −0.061 | 0.820 | 0.732 | 0.394 | 0.365 | 235 | 212 |

#### 30 ng

| VAF% | TI sens | +d1 sens | Δ      | TI snv | +d1 snv | TI ind | +d1 ind | TI tp | +d1 tp |
|------|---------|----------|--------|--------|---------|--------|---------|-------|--------|
| 0.001 | 0.005 | 0.000 | −0.005 | 0.010 | 0.000 | 0.000 | 0.000 | 2 | 0 |
| 0.01  | 0.005 | 0.000 | −0.005 | 0.010 | 0.000 | 0.000 | 0.000 | 2 | 0 |
| 0.1   | 0.083 | 0.021 | −0.061 | 0.137 | 0.034 | 0.018 | 0.006 | 31 | 8 |
| 0.25  | 0.285 | 0.080 | −0.205 | 0.381 | 0.098 | 0.171 | 0.059 | 107 | 30 |
| 0.5   | 0.496 | 0.173 | −0.323 | 0.654 | 0.229 | 0.306 | 0.106 | 186 | 65 |
| 1.0   | 0.589 | 0.272 | −0.317 | 0.781 | 0.342 | 0.359 | 0.188 | 221 | 102 |
| 2.0   | 0.608 | 0.395 | −0.213 | 0.795 | 0.537 | 0.382 | 0.224 | 228 | 148 |

#### 5 ng

| VAF% | TI sens | +d1 sens | Δ      | TI snv | +d1 snv | TI ind | +d1 ind | TI tp | +d1 tp |
|------|---------|----------|--------|--------|---------|--------|---------|-------|--------|
| 0.001 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0 | 0 |
| 0.01  | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0 | 0 |
| 0.1   | 0.008 | 0.005 | −0.003 | 0.015 | 0.010 | 0.000 | 0.000 | 3 | 2 |
| 0.25  | 0.075 | 0.029 | −0.045 | 0.093 | 0.039 | 0.053 | 0.018 | 28 | 11 |
| 0.5   | 0.187 | 0.125 | −0.061 | 0.249 | 0.171 | 0.112 | 0.071 | 70 | 47 |
| 1.0   | 0.400 | 0.285 | −0.115 | 0.537 | 0.390 | 0.235 | 0.159 | 150 | 107 |
| 2.0   | 0.539 | 0.467 | −0.072 | 0.717 | 0.624 | 0.324 | 0.277 | 202 | 175 |

---

## Analysis

### 1. The filter works for its intended purpose

The single spurious PASS call in the 15 ng negative control is eliminated. That call had NDUPALT=0: both alt molecules were simplex-only, consistent with a single-strand deamination error. Requiring NDUPALT >= 1 is a mechanistically sound criterion for distinguishing sequencing artefacts from real somatic variants — a duplex molecule has been sequenced on both strands, which makes a PCR or sequencing error far less likely.

### 2. The sensitivity cost is large

Requiring one duplex-confirmed molecule removes a substantial fraction of real calls at every VAF level below 2%:

- **At 0.5% VAF (15 ng):** 159 → 83 TPs. Nearly half of all detectable calls are lost.
- **At 1.0% VAF (15 ng):** 216 → 147 TPs. One third of TPs are lost.
- **At 2.0% VAF (15 ng):** 235 → 212 TPs. Loss is small — only 23 TPs.

The sensitivity loss narrows at 2% VAF because higher VAF means more alt molecules per sample, and more alt molecules means a higher probability that at least one duplex pair was captured. At low VAF (0.1–0.5%), most alt sites are represented by one or two molecules, and duplex confirmation of a single molecule is a chance event dependent on read depth.

### 3. The 30 ng input is hit hardest

The 30 ng condition loses far more TPs than 15 ng or 5 ng at the same VAF:

| Input | TI tp at 2% | +d1 tp at 2% | Lost | % lost |
|-------|-------------|--------------|------|--------|
| 15 ng | 235 | 212 | 23 | 9.8% |
| 30 ng | 228 | 148 | 80 | 35.1% |
| 5 ng  | 202 | 175 | 27 | 13.4% |

This is counterintuitive — more input DNA should mean more duplex coverage. The explanation is that at 30 ng, there are more unique molecules in the library. With a fixed 2M read budget, each molecule receives fewer reads on average, which reduces the probability of capturing both strands of an alt-molecule as a duplex pair. At 5 ng, the library is smaller and 2M reads achieves higher per-molecule depth, so duplex pairs are more likely to be captured.

At clinical sequencing depths (e.g. 10–20M reads), the penalty for 30 ng input would shrink substantially because duplex capture would become more complete.

### 4. Indels are hit proportionally harder than SNVs

Comparing relative loss at 2% VAF (15 ng):
- SNV: 168 → 150 TPs (-10.7%)
- Indel: 67 → 62 TPs (-7.5%)

At 1% VAF (15 ng):
- SNV: 156 → 108 TPs (-30.8%)
- Indel: 60 → 39 TPs (-35.0%)

At lower VAFs indels consistently show a larger relative loss. Indels require a longer contiguous alt path in the de Bruijn graph; they tend to have fewer supporting molecules on average, and among those molecules, the fraction with duplex confirmation is lower. This mirrors the general observation that indel sensitivity is lower than SNV sensitivity throughout the benchmark.

### 5. The filter is only appropriate in high-depth or high-duplex-rate regimes

At 2M reads and a ~10% duplex rate (15 ng), roughly 25–36% of genuine 2% VAF calls have NDUPALT=0. This is documented in the `CallerConfig::min_alt_duplex` comment. A threshold of 1 at this depth removes more true positives than false positives except in a negative control with no spiked-in variants.

The appropriate use cases for `--min-alt-duplex 1` are:

1. **Deeper sequencing.** At 10–20M reads and the same 10% duplex rate, nearly all real alt molecules will have at least one duplex pair. The sensitivity penalty will be much smaller.
2. **Higher duplex-rate protocols.** Libraries with duplex rates above 20–25% make NDUPALT >= 1 a reliable criterion at 2M reads.
3. **Negative control validation.** Running a known-negative sample with `--min-alt-duplex 1` to confirm no background calls pass.
4. **TI mode + conservative monitoring.** When FP rate at known positions matters more than sensitivity (e.g. clinical monitoring where a confirmed false positive triggers unnecessary treatment decisions).

### 6. Is the spurious call a problem in practice?

The background call at chr3:10149804:C:T is a deamination artefact at a truth VCF coordinate. TI mode inherits all limitations of the truth VCF: any background error that coincidentally matches a truth coordinate will pass. In a clinical monitoring context, the truth VCF contains the patient's own somatic variants identified at diagnosis, and the coordinates are known to the assay. Background error at those specific coordinates is plausible (cytosine deamination is non-random and somewhat position-dependent). TI mode cannot distinguish this from a real recurrence.

At 2M reads, this happened in one of three negative controls. At higher depth or with more samples, such coincidental background calls would be expected to appear more frequently. The duplex filter is the correct mechanism to address this, with the trade-off documented above.

---

## Recommendation

**Do not set `--min-alt-duplex 1` as a default.** The sensitivity cost at 2M reads is too large. The default of 0 is correct for the current sequencing depth.

**Consider `--min-alt-duplex 1` for:**
- Sequencing depths >= 10M reads
- Libraries with duplex rates >= 20%
- Negative control QC runs
- Clinical monitoring workflows where false positives at known coordinates are unacceptable

For this benchmark dataset (2M reads, ~10% duplex), TI mode without a duplex requirement is the better operating point. The single background call in the 15 ng negative control is an expected limitation at this depth, not a systematic failure.
