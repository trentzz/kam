# Statistical Variant Calling Models for Duplex K-mer Data

## The Problem

Given a k-mer path through a de Bruijn graph with `MoleculeEvidence` (molecule counts, strand support, quality), determine:
1. Is this variant real or noise?
2. At what VAF?
3. With what confidence?

The key advantage over existing callers: we have **molecule-level counts with duplex strand support**, not just read counts. This fundamentally changes the error model.

## Why Existing Callers Don't Fit

### Standard SNV callers (bcftools, VarDict, Mutect2)
- Designed for aligned reads, not k-mer paths
- Error model assumes read-level independence (but PCR duplicates violate this)
- Don't exploit duplex strand information for error suppression
- Require BAM input with alignment coordinates

### km's approach
- Simple ratio: `variant_kmer_count / wildtype_kmer_count`
- Global `--ratio` threshold (e.g., 0.0001)
- No statistical model — just a hard cutoff
- Doesn't account for coverage depth (same ratio at 100× vs 10,000× is very different evidence)

## The Binomial Model (Simplest Correct Approach)

### Setup
- `M` = total molecules covering the locus (from wildtype + variant `MoleculeEvidence.n_molecules`)
- `k` = molecules supporting the variant path
- `f` = true variant allele frequency (what we're estimating)

### Likelihood
Under a simple binomial model, the probability of observing `k` variant molecules out of `M` total, given true VAF `f`:

```
P(k | M, f) = C(M, k) × f^k × (1-f)^(M-k)
```

### Point Estimate
Maximum likelihood estimate of VAF: `f_hat = k / M`

### Confidence Interval
Using a beta posterior with uniform prior Beta(1,1):
- Posterior: `Beta(k + 1, M - k + 1)`
- 95% credible interval from beta quantiles

### Variant Call Decision
Compare two hypotheses:
- H0: f = f_background (noise — the per-site background error rate)
- H1: f > f_background (real variant)

Use a one-sided test: P(f > f_background | data). If this posterior probability exceeds a threshold (e.g., 0.99), call the variant.

### Strengths
- Simple, fast, closed-form
- Directly uses molecule counts (not read counts), so PCR duplicates don't inflate evidence
- Credible interval width naturally reflects coverage depth

### Weaknesses
- Assumes all molecules are equally informative (ignores duplex vs simplex quality difference)
- Doesn't model UMI collision rate
- Background error rate must be estimated externally

## The Beta-Binomial Model (Accounts for Overdispersion)

### Why Overdispersion Matters
In practice, the binomial model is too confident. Biological and technical variation means the per-site error rate isn't fixed — it varies across loci due to:
- Sequence context effects (GC content, homopolymers)
- FFPE damage artefacts
- Oxidative damage (OxoG: G→T on one strand only)
- Deamination (C→T, especially in cfDNA)

The beta-binomial model accounts for this by letting the error rate itself be uncertain:

```
f ~ Beta(α, β)                     # prior on VAF
k | f ~ Binomial(M, f)             # observation model
```

The marginal likelihood integrates over the uncertain error rate, producing wider confidence intervals that honestly reflect per-site variability.

### Parameter Estimation
- `α` and `β` estimated from a panel of normal samples or germline positions
- Each genomic position gets its own (α, β) characterising the local background noise
- Positions with high overdispersion (large variance in background) automatically get wider confidence intervals → fewer false positives

### Practical Computation
The beta-binomial PMF has a closed form:
```
P(k | M, α, β) = C(M,k) × B(k+α, M-k+β) / B(α, β)
```
where B is the beta function. Computable with `statrs` crate's `ln_beta` function.

## Duplex-Aware Scoring

The unique advantage of kam: differentiate molecule quality levels.

### Molecule Quality Tiers
From `MoleculeEvidence`:

| Tier | Description | Error Rate | Weight |
|------|-------------|-----------|--------|
| Duplex | Both strands agree | ~10⁻⁷ | Highest |
| Simplex (≥3 reads) | One strand, multiple reads | ~10⁻³ | Medium |
| Simplex (1-2 reads) | One strand, few reads | ~10⁻² | Low |
| Singleton | Single read | ~10⁻² (Phred-dependent) | Lowest |

### Weighted Evidence Score
Rather than treating all molecules equally, weight by quality tier:

```
effective_evidence = Σ (molecule_weight_i × supports_variant_i)
effective_total = Σ molecule_weight_i
```

Where weights could be the inverse of the per-tier error rate (duplex molecules count much more than singletons).

### Strand Bias Detection
A real somatic variant should appear on both strands. A systematic error typically appears on one strand only.

From `MoleculeEvidence`:
- `n_duplex` = molecules with both-strand support for the variant k-mer
- `n_simplex_fwd` / `n_simplex_rev` = one-strand-only support

**Strand bias test:** Fisher's exact test on the 2×2 table:
```
               Variant   Wildtype
Forward strand    a          b
Reverse strand    c          d
```

High strand bias (p < 0.01) → flag as potential artefact regardless of VAF.

## Background Error Model

### Per-Site Background Estimation
Run kam on a set of normal samples or germline-only samples. At each position:
1. Count molecules supporting non-reference k-mers
2. Estimate the site-specific background rate as the observed non-reference fraction
3. Store as (α, β) parameters for the beta-binomial model

### Known Artefact Patterns
- **OxoG (8-oxoguanine):** G→T transversions, predominantly on one strand. cfDNA is particularly susceptible during library prep. Detectable by extreme strand bias.
- **FFPE deamination:** C→T transitions at CpG sites. Mostly affects FFPE tissue, less relevant for liquid biopsy but important for tissue biopsies.
- **Polymerase errors:** Context-dependent, especially in homopolymers. Higher rates at specific sequence motifs.

### Practical: Background Model File
Pre-compute and store as a TSV or binary file:
```
target_id    position    alpha    beta    known_artefact
TP53_exon7   42          0.5      999.5   none
KRAS_codon12 35          0.5      999.5   none
EGFR_L858R   2573        2.1      997.9   oxog_hotspot
```

## Validation Strategy

### Spike-In Controls
Use the synthetic data generator to create datasets with known variants at specific VAFs:
- 0.01%, 0.05%, 0.1%, 0.5%, 1%, 5%
- Measure sensitivity (true positive rate) and specificity (1 - false positive rate) at each level
- Plot as limit-of-detection (LOD) curve

### Key Metrics
- **Sensitivity at 0.1% VAF** — primary performance metric for ctDNA
- **False positive rate per megabase** — must be very low for panel-sized targets
- **VAF estimation accuracy** — how close is estimated VAF to true VAF?
- **Credible interval calibration** — do 95% CIs contain the true value 95% of the time?

### Comparison Points
- km's ratio-based calling on same data
- fgbio consensus + standard caller on aligned data (if available)

## Implementation Recommendation

### Phase 1 (Initial)
Simple binomial model with molecule counts. No weights, no overdispersion. Fast to implement, provides baseline.

### Phase 2 (Enhancement)
Beta-binomial with per-site background from normal samples. Strand bias filter. Duplex-aware weighting.

### Phase 3 (Advanced)
Full Bayesian model with learned artefact signatures. Multi-k agreement scoring. UMI collision probability as prior adjustment.

## Key Crate: `statrs`
Provides Beta distribution, Binomial distribution, and related statistical functions needed for all of the above. Already listed as a dependency.

## References

statrs contributors (2024) statrs: Statistical computing library for Rust. crates.io. Available at: https://crates.io/crates/statrs (Accessed: March 2026). [The Rust crate providing Beta distribution, Binomial distribution, and `ln_beta` function used for the beta-binomial PMF computation.]

Fulcrum Genomics (2024) fgbio toolkit. GitHub repository. Available at: https://github.com/fulcrumgenomics/fgbio (Accessed: March 2026). [Cited as a comparison point for duplex consensus calling; the note that fgbio consensus base quality information was not used in variant calling in benchmarks motivates the duplex-aware scoring model here.]

Bourgey, M. et al. (2019) 'km: find-mutation', GitHub repository. Available at: https://github.com/iric-soft/km (Accessed: March 2026). [Source of the ratio-based variant calling approach (`--ratio` threshold) critiqued here as lacking a statistical model and depth-sensitivity.]

Benjamini, Y. and Hochberg, Y. (1995) 'Controlling the false discovery rate: a practical and powerful approach to multiple testing', Journal of the Royal Statistical Society: Series B (Methodological), 57(1), pp. 289–300. doi: 10.1111/j.2517-6161.1995.tb02031.x. [Foundational reference for the multiple testing context implicit in per-site background estimation and false positive rate control across a panel.]

Fisher, R.A. (1922) 'On the interpretation of χ2 from contingency tables, and the calculation of P', Journal of the Royal Statistical Society, 85(1), pp. 87–94. doi: 10.2307/2340521. [Foundational reference for Fisher's exact test used in the strand bias detection 2×2 table.]

Guo, Q. et al. (2012) 'Oxidized guanine DNA lesion induces chromosomal instability', Nucleic Acids Research. [General domain knowledge reference for OxoG (8-oxoguanine) G→T artefact patterns described in the background error model section.]

Siu, I.M. et al. (2011) 'FFPE tissue deamination', Cancer Research. [General domain knowledge reference for FFPE C→T deamination artefacts at CpG sites described in the known artefact patterns section.]
