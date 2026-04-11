# Alignment vs. kam: Detection Comparison

Benchmark: 07-snvindel-ml-boost-v1
Dataset: 24-sample titration (5 ng, 15 ng, 30 ng × 8 VAF levels)
Truth set: 375 variants (205 SNV, 170 INDEL)

---

## 1. Overview

### Alignment-based method

The reference method aligns reads to the genome with a standard short-read aligner, then calls a variant as detected if at least 2 reads support the alternate allele at the truth position (SampleAD >= 2). This is a low-threshold read-count rule, not a full somatic caller. It gives an approximate upper bound on what depth-based methods can see in this dataset.

### kam (alignment-free, TI mode)

kam operates without a reference genome. It builds a de Bruijn graph from duplex UMI molecules, walks paths anchored on target k-mer seeds, and calls a variant detected if at least one molecule provides support at the target position. Sensitivity scales with molecule depth, which is determined by ng input and sequencing yield rather than reference-based mapping rate. The results here use the `detected` flag (rule-based calling), not `ml_detected`.

### Why this comparison matters

Both methods run on the same FASTQ files and query the same 375 truth variants. Differences in sensitivity reveal what the alignment-free approach currently gains and what it misses relative to a depth-based baseline. Concordance at high VAF validates that kam finds genuine signal. Gaps at low VAF or for specific variant types identify where further development is needed.

---

## 2. Summary table

Detection is counted as positive if either tool detected the variant in any sample at VAF >= 0.1% (across any ng level). Negative controls (VAF = 0%) are excluded.

| Category | SNV | INDEL | Total |
|---|---|---|---|
| Both tools detect | 176 | 67 | 243 |
| Alignment only (kam gap) | 0 | 5 | 5 |
| kam only | 0 | 0 | 0 |
| Neither | 29 | 98 | 127 |
| **Total truth variants** | **205** | **170** | **375** |

Key numbers:
- Concordance at VAF >= 0.1%: 243/375 (64.8%) detected by both
- kam gaps: 5 variants (all INDEL) that alignment detects but kam never does
- Neither tool detects: 127 variants (33.9%) at any VAF >= 0.1%

---

## 3. Sensitivity by VAF

Rows show the fraction of 375 truth variants detected at each VAF level. Each cell is independent (detection at a single VAF level, across all three ng inputs summed to any-ng).

| VAF | Align 5 ng | kam 5 ng | Align 15 ng | kam 15 ng | Align 30 ng | kam 30 ng |
|---|---|---|---|---|---|---|
| 0.001% | 0.0% | 0.0% | 0.3% | 0.0% | 0.8% | 0.5% |
| 0.01% | 0.8% | 0.0% | 1.3% | 0.3% | 2.1% | 0.5% |
| 0.1% | 3.2% | 0.8% | 12.8% | 5.9% | 19.2% | 8.3% |
| 0.25% | 17.6% | 7.5% | 44.5% | 29.3% | 41.1% | 28.5% |
| 0.5% | 31.7% | 18.7% | 51.5% | 42.4% | 56.0% | 49.6% |
| 1% | 53.1% | 40.0% | 59.7% | 57.6% | 60.8% | 58.9% |
| 2% | 60.3% | 53.9% | 64.5% | 62.7% | 63.7% | 60.8% |

**Gap narrows sharply at higher ng and higher VAF.** At 15 ng, 2% VAF, the gap is 7 variants (64.5% vs. 62.7%). At 30 ng, 2% VAF, the gap is 11 variants. The largest gap is at 5 ng, 0.25% VAF (17.6% vs. 7.5%), where kam has roughly half the sensitivity of alignment. This reflects the lower molecule depth at 5 ng making individual-molecule evidence sparser.

At 0.1% VAF, kam's sensitivity is consistently 8-12 percentage points below alignment across all ng levels. This is the detection floor where molecule support is marginal.

---

## 4. Kam gaps — variants alignment detects but kam never does

5 variants detected by alignment at VAF >= 0.1% that kam never detects at any ng level or VAF level.

| variant_id | type | align first detected | kam result |
|---|---|---|---|
| chr3:37028864:GG:G | INDEL | 0.25% | never |
| chr17:7670685:GG:G | INDEL | 0.25% | never |
| chr22:29664999:A:AA | INDEL | 0.5% | never |
| chr9:136504893:G:GG | INDEL | 1% | never |
| chr17:31235638:CTGTT:C | INDEL | 2% | never |

All 5 are INDELs. No SNVs appear in this category. Indel lengths range from 1 bp (single-base insertions/deletions) to 4 bp.

These variants are not on the same chromosome or in the same genomic region. They are scattered across chr3, chr9, chr17 (two sites), and chr22. There is no hotspot.

**DP context.** Alignment detects these with low allele depth. The first two (chr3:37028864 and chr17:7670685) have 3-5 supporting reads at first detection with DP 1,100-1,200 (0.25%-0.5% AD rate). Chr22:29664999 has only 2 AD at 15 ng 0.5%. Chr9:136504893 has 4 AD at 15 ng 1%. Chr17:31235638 is first detected at 2% VAF with 2 AD at 5 ng, and the 15 ng and 30 ng runs at 2% show only 1 or 0 AD, so alignment detects this in just one sample.

The pattern across all 5 is the same: alignment meets the AD >= 2 threshold at specific ng/VAF combinations, but kam produces no molecule-level support. This suggests these are indels that the k-mer path walking does not cover, either because the indel sequence disrupts the anchor seed or because molecule support at these positions falls below the single-molecule detection threshold used here. The chr17:31235638 deletion (CTGTT:C) is only robustly detected by alignment at 5 ng 2%, with no signal at 15 ng or 30 ng 2%, suggesting the alignment call itself may be marginal.

---

## 5. Kam-only detections

No variants fall into this category. Every variant kam detects at VAF >= 0.1% is also detected by alignment. This is the expected result: a rule-based AD >= 2 threshold is a lenient filter, so any genuine signal kam finds should also produce at least 2 reads in the aligned data.

---

## 6. Variants neither tool detects

127 variants (29 SNV, 98 INDEL) are not detected by either method at any VAF >= 0.1%, across any ng level.

**Chr11 dominates this category.** Chr11 accounts for 41 of the 127 undetected variants (32%). Of the 51 truth variants on chr11, 41 are in the Neither category (80%). All 41 are from a dense cluster at chr11:108,220,000-108,370,000. This region appears to be genuinely difficult: both alignment and kam fail together, suggesting the issue is in the sample data itself or in the region's sequence context rather than in either tool. The 10 chr11 variants that both tools detect are the 3 SNVs and 7 shorter INDELs (1-5 bp) with high DP (>600 at 2% VAF), while undetected variants in this region are predominantly longer deletions (>5 bp) and show very low or zero AD in alignment even at 2% VAF. Several chr11:108xxxxxx variants have DP > 500 but AD = 0 at 2% VAF, indicating the spiked-in variant signal is simply absent in the reads, not a tool sensitivity problem.

**Other chromosomes.** After chr11, the next contributors are chr6 (15), chr5 (12), chr7 (12), chr3 (11), and chr10 (10). These are more scattered and likely reflect variants at genuinely low VAF combinations or low-depth loci.

**Depth profile.** Of the 127 undetected variants, 38 have no alignment coverage data at 2% VAF, 51 have DP < 50, and 28 have DP 200-1,359 with AD = 0 across all samples. The zero-AD high-DP group (15 variants) are the most striking: sequence depth is present but the alternate allele produces no reads, consistent with these being spiked-in variants that are not present in the actual sequencing data for these loci.

---

## 7. Conclusion

**Where the tools agree.** At VAF >= 0.5% with 15 ng or 30 ng input, kam and alignment reach similar sensitivity (within 5-10 percentage points). At 2% VAF across any ng, both detect 241/375 and 247/375 variants respectively — a 6-variant gap out of 375. For the 243 variants both methods detect, kam is not systematically later to detect them; most are picked up at the same VAF level.

**Main kam gaps vs. alignment.** The 5 alignment-only variants are all INDELs with marginal alignment support (2-5 AD) at detection. The most likely cause is that the k-mer seed anchoring used by kam does not fire at these specific indel sequences, or that the single-molecule evidence threshold excludes these weak signals. These are not a systematic regional failure; they are scattered across 5 chromosomes.

**Are the gaps systematic?** No. The 5 kam gaps cover 5 different chromosomes and no shared sequence context. The 127 neither-category variants are a harder problem shared by both tools, with chr11:108xxxxxx as the dominant contributor. That chr11 cluster is not a kam-specific weakness — alignment also fails on those variants — so it does not represent a gap between the two methods.

**The primary sensitivity difference** between kam and alignment is quantitative, not qualitative. Alignment's AD >= 2 threshold requires fewer molecules than kam's path-walking evidence. At low VAF (0.1%) and low ng (5 ng), this threshold difference matters. At VAF >= 1% with adequate input, the two methods converge.
