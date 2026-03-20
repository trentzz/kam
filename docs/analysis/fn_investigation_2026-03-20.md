# False Negative Investigation — 30ng 2% VAF, 2M reads

**Date:** 2026-03-20
**Sample:** TWIST_STDV2_30ng_VAF_2pc (2,000,000 read pairs)
**Molecules assembled:** 397,613 (26,230 duplex, 6.6%)
**Pipeline run:** kam v0.1, all three fixes applied (strand bias min, duplex orientation, mean_duplex filter)

---

## Summary

375 truth variants were evaluated (205 SNPs, 170 indels).
224 were detected as PASS (59.7% sensitivity).
151 were missed (false negatives).

| Category | Count | SNP | INDEL |
|---|---|---|---|
| TP (PASS) | 224 | 160 | 64 |
| FN: no alt path in graph | 112 | 23 | 89 |
| FN: path exists but variant not found | 27 | 14 | 13 |
| FN: LowConfidence (n_alt = 1) | 10 | 7 | 3 |
| FN: StrandBias | 2 | 2 | 0 |
| FN: LowDuplex | 0 | 0 | 0 |

The dominant cause of missed variants is the **no-path category (112 out of 151 FN)**, followed by **path-exists-but-not-found (27)** and **single-molecule detections blocked by the n_alt ≥ 2 threshold (10)**.

---

## Cause 1: No alt path for target (112 FN)

375 targets were processed. 263 targets produced at least one alt path (268 targets with variants per the log — slight discrepancy because 5 targets produced paths that did not score as truth variant matches). 112 targets produced zero alt paths.

All 112 FN in this category come from targets with no alt paths. **89 of the 112 are indels (79%).**

### Why no alt path?

The de Bruijn graph is built from k-mers observed in the assembled molecules. For a variant path to be walkable:

1. The alt-specific k-mers must be present in the graph (requires ≥ `min_molecules` coverage of the alt k-mers).
2. The end-anchor k-mer must be reachable from the alt region.

For SNPs (23 cases): the most likely cause is **stochastic molecule dropout**. At 2% VAF with 397K molecules total and 375 targets, the expected alt molecule count per target is approximately `397,613 × 0.02 / 375 ≈ 21` molecules. After UMI grouping, collisions, and k-mer filtering, the effective count reaching the graph is lower. Some targets simply had no alt molecules survive the pipeline.

For indels (89 cases): an additional structural problem compounds dropout. The end-anchor k-mer of the reference path **shifts by the indel size** relative to the genomic sequence. A deletion of length ℓ moves the end anchor ℓ bases earlier; an insertion moves it ℓ bases later. The shifted anchor k-mer must be:
- present in the sequencing reads from alt-supporting molecules, AND
- covered by enough molecules to pass the `min_molecules = 1` threshold in the graph.

At 2% VAF, even when alt molecules are present, their reads may not span the shifted anchor position — particularly for large indels where the shift can be 10–30 bases. This is why large indels are almost never detected (see Cause 4).

**Ref molecule coverage** for the 263 targets that did produce paths ranged from 1 to 1,189, with a median of 539 molecules. Zero-path targets are likely at the lower end of this distribution, though we cannot confirm without per-target molecule tracking.

---

## Cause 2: Path exists but truth variant not found (27 FN)

27 truth variants are in targets that produced at least one alt path, but the specific truth variant was not among the called variants (PASS or filtered). These break down into 14 SNP and 13 INDEL.

### Likely sub-causes

**a. Insufficient alt k-mer coverage despite partial graph presence.** A target might have one well-covered variant (e.g., a common germline SNP) that produces a walkable path, while the rare somatic variant at 2% VAF does not have enough k-mers in the graph to form a complete path. The DFS succeeds for the higher-coverage path but never reaches the rare variant.

**b. Path combinatorial explosion or max-paths limit.** Some targets in this category have 3–4 other calls (including PASS calls), suggesting complex graphs. The DFS has a max-paths limit (`max_paths = 8` by default). If the DFS fills the limit with higher-evidence paths before encountering the truth variant's path, the truth variant is never scored.

**c. Normalisation mismatch.** For indels, the left-normalisation applied by the scoring script and by kam must agree exactly. If kam calls the variant at a slightly different position due to a normalisation edge case, the truth-matching fails even though the variant was found.

Examples of interest from this category:

| Variant | Type | Target calls | PASS calls | Note |
|---|---|---|---|---|
| chr3:37028864 GG>G | INDEL 1bp | 6 | 3 | Multi-variant target |
| chr11:108284389 CAGAGACA>C | INDEL 7bp | 3 | 0 | Complex deletion |
| chr11:108301670 GTTACCTGT>G | INDEL 8bp | 2 | 0 | Large deletion |
| chr11:108317458 TAGAAGAA>TAGTC | INDEL 3bp | 1 | 0 | Complex MNV-indel |
| chr22:29664999 A>AA | INDEL 1bp | 2 | 1 | Another variant PASS'd |

The chr11 cluster (5 variants across multiple chr11 targets) suggests that targets in certain high-indel gene regions (possibly MLL2/KMT2D) are consistently difficult.

---

## Cause 3: LowConfidence (n_alt = 1) — 10 FN

10 truth variants were detected as a path with exactly 1 supporting molecule but were blocked by the `min_alt_molecules = 2` threshold. The posterior confidence filter (≥ 0.99) is not the bottleneck here — confidence ranges 0.87–1.00 for these 10, showing the molecule count is the binding constraint.

| Variant | Type | n_alt | n_ref | VAF | Confidence |
|---|---|---|---|---|---|
| chr3:138946321 G>C | SNP | 1 | 62 | 0.0159 | 0.983 |
| chr3:179218303 G>A | SNP | 1 | 586 | 0.0017 | 0.869 |
| chr4:54280374 A>G | SNP | 1 | 351 | 0.0028 | 0.916 |
| chr6:117311094 C>A | SNP | 1 | 33 | 0.0294 | 0.991 |
| chr7:55019353 G>A | SNP | 1 | 13 | 0.0714 | 0.996 |
| chr9:77794572 T>G | SNP | 1 | 392 | 0.0025 | 0.907 |
| chr10:108256298 T>TAAAAA | INS | 1 | 570 | 0.0018 | 0.872 |
| chr11:108331534 G>GA | INS | 1 | 403 | 0.0025 | 0.905 |
| chr11:108365362 GTCT>G | DEL | 1 | 393 | 0.0025 | 0.907 |
| chr19:10491905 C>A | SNP | 1 | 220 | 0.0045 | 0.945 |

The chr7:55019353 case (EGFR T790M) is notable: observed VAF of 7.1% with only 13 ref molecules — this target has very low coverage. The single alt molecule has 99.6% confidence.

The three chr11 cases (108256298, 108331534, 108365362) are a cluster in the MLL2/KMT2D gene with low per-molecule coverage.

**Recovery options:**
- Lower `min_alt_molecules` to 1 when `n_duplex_alt ≥ 1`. Duplex confirmation on a single molecule is strong evidence.
- At current defaults, lowering to 1 unconditionally would recover these 10 at the cost of more single-molecule false positives.

---

## Cause 4: StrandBias — 2 FN

Two SNPs were detected with strong molecule support but flagged as StrandBias:

| Variant | Gene | n_alt | n_ref | n_duplex (mean) | sb_p | VAF |
|---|---|---|---|---|---|---|
| chr12:25245347 C>T | KRAS | 11 | 860 | 46 | 0.0064 | 0.0126 |
| chr22:21772875 C>T | CRKL | 17 | 715 | 39 | 0.0050 | 0.0232 |

Both have n_simplex_alt = 0 in the output (derived as `max(0, n_alt − mean_duplex)`), but the Fisher test still gives p < 0.01. This is not a contradiction.

### Root cause

The strand bias test uses `min_simplex_fwd` and `min_simplex_rev` from `PathEvidence` — the **minimum per-direction simplex count across all k-mers in the alt path**. For variant-specific k-mers (the positions that differ from reference), the only contributing molecules are the alt-supporting molecules. If those alt molecules happen to be all forward-strand simplex (no reverse-strand coverage), then `min_simplex_fwd = n_alt` and `min_simplex_rev = 0` at those k-mers, giving the Fisher test input (11, 0, ref_fwd, ref_rev) or (17, 0, ref_fwd, ref_rev).

The Fisher test correctly identifies this as strand-imbalanced. But at n_alt = 11–17, the probability of observing all alt molecules from one strand by pure sampling chance (binomial(11, 0.5) ≤ 1 success) is approximately 0.05–0.005%, which falls just below the p < 0.01 threshold. These are **type I errors driven by small sample size**, not genuine strand bias.

The mean_duplex displayed (39–46) comes from anchor k-mers that are shared with the reference path and are covered by many duplex reference molecules, artificially inflating the displayed duplex count.

**Recommendation:** Skip or relax the strand bias test when `n_alt < 20`. At fewer than 20 alt molecules, the power to detect genuine strand bias is limited and sampling artefacts are common. Alternatively, apply a Bonferroni correction for low-count variants.

---

## Cause 5: LowDuplex — 0 FN (filter is ineffective)

Zero truth variants were blocked by the LowDuplex filter. Investigation shows that 0 calls in the entire run were labelled LowDuplex.

**Why:** The `n_duplex_alt` used in the filter is `alt_evidence.mean_duplex.round() as u32`. Because the alt path contains ~70 k-mers including anchor k-mers shared with reference, and those anchor k-mers are covered by many duplex reference molecules, `mean_duplex` is nearly always ≥ 1 even for spurious low-evidence calls. The LowDuplex filter therefore never triggers.

**Implication:** The min_alt_duplex = 1 constraint provides no additional discrimination over the min_alt_molecules = 2 constraint, because every alt path that clears the molecule threshold will also have mean_duplex ≥ 1.

**Fix:** Change `n_duplex_alt` from `mean_duplex` to the minimum duplex count across variant-specific k-mers only (excluding anchor k-mers), or use a different metric such as `min_duplex` across the central k-mers of the path.

---

## Indel size effect

TP and FN indels show a clear size cut-off:

| Indel size (bp) | TP | FN |
|---|---|---|
| 1 | 29 | 15 |
| 2 | 11 | 9 |
| 3 | 11 | 5 |
| 4 | 3 | 10 |
| 5 | 1 | 7 |
| 6 | 3 | 4 |
| 7 | 4 | 4 |
| 8 | 2 | 4 |
| > 8 | 0 | 45 |

No indel larger than 8 bp is ever detected. For indels of 1–3 bp, roughly half are still missed. Both observations are consistent with the anchor shift problem: larger indels move the end anchor further from the reference position, requiring reads to cover a wider shifted region. At 2% VAF, variant-molecule coverage of the shifted anchor position decreases with indel size.

For 1-bp indels, 15/44 = 34% are missed despite the smallest possible anchor shift. This points to a secondary cause: the 1-bp indel k-mer paths are in a repeat or low-complexity context where the shifted anchor k-mer is non-unique (present in the reference path and elsewhere), causing path enumeration to miss or confuse the alt path.

---

## Prioritised recommendations

| Priority | Fix | FN recovered | Complexity |
|---|---|---|---|
| 1 | Increase depth (hash-partition UMI grouping) | ~30–50 at 5M reads | High |
| 2 | Lower min_alt_molecules to 1 when min_duplex ≥ 1 | 10 | Low |
| 3 | Skip strand bias test when n_alt < 20 | 2 | Low |
| 4 | Fix min_duplex computation (use variant-specific k-mers only) | Unclear | Medium |
| 5 | Investigate path-not-found category (27 FN) more deeply | Up to 27 | Medium |

The biggest single win is depth: the 112 zero-path FN represent variants where the alt molecules simply weren't assembled or their k-mers weren't covered. More reads means more alt molecules per target. The O(n²) UMI grouping step is the current bottleneck limiting practical depth.

---

## Data

- Sample: `TWIST_STDV2_30ng_VAF_2pc_22KVL2LT3`
- Output: `/tmp/kam_investigate/out/variants.tsv` (1,216 variants: 316 PASS, 889 LowConfidence, 11 StrandBias)
- Truth: `benchmarking/scripts/truth_variants.vcf` (205 SNP, 170 INDEL)
