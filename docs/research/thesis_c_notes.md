# Thesis C Notes: k-mer Variant Detection for Duplex UMI Sequencing

**Source**: Trent Zeng, UNSW, November 2025
**Title**: "Implementing a k-mer-based variant detection method for UMI-based duplex sequencing"
**Location**: `/home/tzeng/tfiles/unsw/thesis/Thesis C/Thesis_C_Final_Report_Masked.pdf`

---

## Overview

The thesis evaluated an alignment-free variant detection pipeline (HUMID + Jellyfish + km) on the same Twist cfDNA Pan-Cancer Reference Standard v2 (STDV2) titration series used for kam benchmarking. The key difference from kam: the thesis pipeline used **tumour-informed filtering**, which is why it achieved near-zero false positives.

---

## Pipeline

```
FASTQ → HUMID (deduplication) → Jellyfish (k-mer counting) → km (path walking) → tumour-informed filter → calls
```

Chapter 3 compared three deduplication strategies on the same samples:

| Strategy   | Description                              |
|------------|------------------------------------------|
| Non-dedup  | No deduplication, raw reads              |
| REDUX      | Alignment-based deduplication            |
| HUMID      | Alignment-free UMI-based deduplication   |

---

## Dataset

Same Twist STDV2 titration series: 5 ng, 15 ng, 30 ng input DNA, VAF 0%–2%.

**Truth set differs from kam**: thesis uses 456 variants (225 SNVs, 49 insertions, 167 deletions, 15 SVs); kam uses 375 (205 SNVs, 170 indels). Both are from the same samples but with different truth set versions.

---

## Tumour-Informed Filtering (Section 3.4)

This is the mechanism behind near-zero false positives.

**How it works**: For each known target mutation (chr, pos, ref, alt), generate both the reference and expected alt k-mer sequences. Run km. Keep only calls where the called alt sequence matches the **expected alt allele exactly**. Discard everything else.

**Quality thresholds applied** (Table 3.3):

| Filter       | Threshold |
|-------------|-----------|
| rVAF        | > 0       |
| expression  | > 0       |
| min_coverage| ≥ 3 (all k-mers on alt path) |

**Why this gives near-zero FPs**: Background cfDNA biological variants, germline SNPs, and CHIP events are invisible. They would only appear as FPs if they landed at the exact same (chrom, pos, ref, alt) as a target mutation, which essentially never happens.

**Critical distinction**: This is **monitoring mode**, not **discovery mode**. The caller only confirms whether a pre-specified known variant is present. It cannot find novel somatic variants.

---

## Results at 2% VAF (30 ng)

### Confusion matrix (Figure 3.8a — HUMID, 30 ng, 2% VAF)

The axes compare km calls vs alignment-based calls:
- 151 TP (km found it, alignment found it)
- 0 extra km calls not found by alignment ("FP" in this sense)
- 18 FN (alignment found it, km did not)

Non-dedup at same conditions: 30 TP, 139 FP — deduplication is essential.

### Sensitivity

| Strategy  | SNV sensitivity | INDEL sensitivity |
|-----------|----------------|-------------------|
| HUMID     | 77%            | 38%               |
| REDUX     | 78%            | 38%               |
| Non-dedup | 39%            | 36%               |

Observation: HUMID and REDUX perform equivalently for both SNVs and INDELs. Non-dedup collapses SNV sensitivity but barely affects INDEL sensitivity — indel calling depends on structural path coverage, not read-level error rates.

---

## Chapter 4: Clinical Validation

Applied the same HUMID+km pipeline to leukemia patient personalised panels (AML patients, deep sequencing).

- 74 SNV queries across 10 patients
- **73.7% overall sensitivity**
- Context: each query was a known somatic variant from patient tumour profiling — pure monitoring mode

---

## Section 5.2.1: KAM Workflow

The thesis described a KAM (K-mer Analysis Modules) workflow as a planned next step. This is the conceptual precursor to the current kam Rust implementation. The thesis was written before the current kam was built.

---

## Implications for kam

### Why kam has more false positives than the thesis

kam operates in **discovery mode**: it finds all alternative paths in every 100 bp target window and reports them. Background biological cfDNA variants (germline SNPs, CHIP, somatic variants in normal tissue) are real biological signal and pass all quality filters.

The thesis avoided these by design — it never looked for them.

### Path to near-zero FPs in kam

Add a `--target-variants VCF` flag implementing **monitoring mode**:

1. Accept a VCF of known expected variants (chrom, pos, ref, alt)
2. After variant calling, filter to only calls where (chrom, pos, ref, alt) matches an entry in the VCF
3. Discard all other discovered variants

This converts kam from a somatic discovery tool to a confirmatory monitoring tool, replicating the tumour-informed filtering behaviour. Expected outcome: FP count drops to near zero (only variants coincidentally matching a target allele at the same position).

### min_coverage≥3 equivalent

The thesis required every k-mer on the alt path to appear ≥3 times after deduplication. kam currently uses `n_alt ≥ 2` molecules. The thesis threshold corresponds roughly to `--min-alt-molecules 3` in kam. This is a secondary filter; tumour-informed filtering is the primary mechanism.

---

## Key Differences: Thesis vs kam

| Aspect               | Thesis (HUMID+km)              | kam                            |
|----------------------|-------------------------------|-------------------------------|
| Mode                 | Monitoring (known variants)    | Discovery (all variants)      |
| FP source            | Near-zero (allele-exact match) | Background biological signal  |
| FP count at 2% VAF   | ~0                             | 35–72 per sample (floor)      |
| Truth set size       | 456 variants                  | 375 variants                  |
| SNV sensitivity      | 77% (30 ng, 2% VAF)           | 69–80% (2 M reads, 2% VAF)    |
| INDEL sensitivity    | 38% (30 ng, 2% VAF)           | 32–39% (2 M reads, 2% VAF)    |
| Calling unit         | k-mer counts (read-level)     | Molecule counts               |
| Duplex evidence      | Lost after HUMID              | Preserved throughout          |
| Runtime              | Multi-tool chain              | 19–35 s per sample, single core|
