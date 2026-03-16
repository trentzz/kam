# Consensus Calling Algorithms for Duplex UMI Families

## What Consensus Calling Is

When multiple reads come from the same original DNA molecule (identified by their shared UMI), they should all carry the same sequence — any differences between them are sequencing or PCR errors. Consensus calling collapses a family of reads into a single "consensus" sequence that represents the best guess at the true original molecule sequence, along with a per-base confidence estimate.

This is the core error-correction step in duplex sequencing. Without it, you're limited by the per-read error rate (~10⁻² to 10⁻³). With single-strand consensus, you get ~10⁻³ to 10⁻⁴. With duplex consensus (both strands agree), you get ~10⁻⁷. This is what makes 0.1% VAF detection possible.

## The Three Levels of Consensus

### Level 1: Single-Strand Consensus (SSC)

Take all reads from one strand of a molecule family (e.g., all forward-strand reads) and collapse them.

```
Read 1 (Q30): A C G T A C G T
Read 2 (Q25): A C G T A C G T
Read 3 (Q20): A C G T A C A T    ← error at position 7
Read 4 (Q28): A C G T A C G T
                              ↓
SSC:          A C G T A C G T    (position 7: 3 votes G vs 1 vote A → G wins)
```

### Level 2: Duplex Consensus (DC)

Take the forward-strand SSC and the reverse-strand SSC and compare them. Positions where both agree have extremely high confidence. Positions where they disagree are ambiguous.

```
SSC forward:  A C G T A C G T
SSC reverse:  A C G T A C G T    ← both strands agree everywhere
                              ↓
Duplex:       A C G T A C G T    (error probability ~10⁻⁷ at each position)
```

### Level 3: Handling Disagreements

When SSC forward and SSC reverse disagree at a position:

```
SSC forward:  A C G T A C G T
SSC reverse:  A C G T A C A T    ← disagree at position 7
```

Options:
- **Strict:** mark position as N (unknown) — safest, loses information
- **Quality-weighted:** pick the base from whichever SSC had higher confidence at that position
- **Flag:** keep both and let downstream decide — report as a "disputed" position

## Algorithm 1: Simple Majority Vote

At each position, the base that appears most frequently across the family wins.

```
Position 7: G G G A G → G wins (4 vs 1)
```

**Pros:** Simple, fast, easy to understand.
**Cons:** Ignores quality. A single Q40 read and four Q5 reads would still pick the Q5 consensus. At small family sizes (2-3 reads), a single error can flip the consensus.

## Algorithm 2: Quality-Weighted Majority Vote — RECOMMENDED STARTING POINT

At each position, each read's base vote is weighted by its Phred quality score (converted to error probability). The base with the highest total weight wins.

### The Math

For each base b at position i across reads in the family:

```
weight(b, i) = Σ (1 - error_prob(read_j, i))    for all reads j where base_j[i] == b
```

Where `error_prob = 10^(-Q/10)` converts Phred to probability.

The consensus base is `argmax_b weight(b, i)`.

### Example

Position 7, family of 4 reads:
```
Read 1: G at Q30 → weight = 1 - 0.001 = 0.999
Read 2: G at Q25 → weight = 1 - 0.003 = 0.997
Read 3: A at Q10 → weight = 1 - 0.100 = 0.900
Read 4: G at Q28 → weight = 1 - 0.0016 = 0.998

G total weight: 0.999 + 0.997 + 0.998 = 2.994
A total weight: 0.900

Consensus: G (overwhelmingly)
```

Now consider a trickier case:
```
Read 1: G at Q10 → weight = 0.900
Read 2: A at Q35 → weight = 0.9997
Read 3: G at Q8  → weight = 0.842

G total weight: 0.900 + 0.842 = 1.742
A total weight: 0.9997

Consensus: G (majority wins despite lower quality)
```

This is where simple majority vote and quality-weighted vote can give different answers. With only 2 reads at Q8-Q10 vs 1 read at Q35, the majority vote picks G but the quality evidence favours A. The quality-weighted vote still picks G here because the weights are additive — but if the Q35 read were the only high-quality read, it might be right.

### Consensus Quality Estimation

After picking the consensus base, estimate how confident we are:

```
p_correct = weight(consensus_base) / Σ weight(all_bases)
consensus_error_prob = 1 - p_correct
consensus_phred = -10 × log10(consensus_error_prob)
```

For the first example: `p_correct = 2.994 / 3.894 = 0.769` — but this underestimates confidence because it treats the weights as probabilities, which they aren't quite. A more principled approach uses the Bayesian model.

**Pros:** Accounts for base quality, fast to compute, handles the common case well.
**Cons:** Additive weighting can be dominated by many low-quality reads outvoting fewer high-quality reads. Quality estimation is approximate.

## Algorithm 3: Bayesian Posterior

Treat consensus calling as Bayesian inference: given the observed bases and their qualities, what is the posterior probability that the true base is A, C, G, or T?

### The Math

Prior: uniform over {A, C, G, T} → P(true = b) = 0.25 for each b.

Likelihood for each read j at position i, given true base b:
```
P(observed_j | true = b) = {
    1 - error_prob_j    if observed_j == b    (correct read)
    error_prob_j / 3    if observed_j != b    (error to one of 3 other bases)
}
```

Posterior (by Bayes' theorem):
```
P(true = b | all reads) ∝ P(true = b) × Π_j P(observed_j | true = b)
```

The consensus base is `argmax_b P(true = b | all reads)`.

### Example

Same position 7, family of 4 reads:
```
Read 1: G at Q30 (error = 0.001)
Read 2: G at Q25 (error = 0.003)
Read 3: A at Q10 (error = 0.100)
Read 4: G at Q28 (error = 0.0016)

P(true=G | data) ∝ 0.25 × (0.999) × (0.997) × (0.100/3) × (0.998)
                 = 0.25 × 0.999 × 0.997 × 0.0333 × 0.998
                 = 0.00829

P(true=A | data) ∝ 0.25 × (0.001/3) × (0.003/3) × (0.900) × (0.0016/3)
                 = 0.25 × 0.000333 × 0.001 × 0.900 × 0.000533
                 = 0.0000000399

Normalise: P(G) ≈ 0.99999+, P(A) ≈ 0.000005
```

The Bayesian approach gives a much sharper answer because it multiplies probabilities (each supporting read contributes multiplicatively, not additively).

### Consensus Quality from Posterior

```
consensus_error_prob = 1 - P(consensus_base | data)
consensus_phred = -10 × log10(consensus_error_prob)
```

This naturally produces very high Phred scores for well-supported consensus bases and properly low scores for ambiguous positions.

**Pros:** Mathematically principled. Produces calibrated confidence estimates. Handles edge cases (mixed quality, small families) correctly. Naturally multiplicative — each agreeing read dramatically increases confidence.
**Cons:** Slightly more computation (log-space multiplications). Assumes independence between reads (violated if errors are systematic, e.g., same position always miscalled due to sequence context). The uniform prior may not be optimal.

### Practical Note: Log-Space Computation

To avoid floating-point underflow with many reads, compute in log space:

```rust
let log_posterior: [f64; 4] = [0.0; 4]; // A, C, G, T
for read in family {
    let error_prob = phred_to_prob(read.quality[pos]);
    for base in 0..4 {
        if read.sequence[pos] == BASES[base] {
            log_posterior[base] += (1.0 - error_prob).ln();
        } else {
            log_posterior[base] += (error_prob / 3.0).ln();
        }
    }
}
// Find argmax, then compute posterior probabilities via log-sum-exp
```

## Algorithm 4: fgbio's Approach

For comparison, fgbio's `CallDuplexConsensusReads` uses a quality-weighted approach similar to Algorithm 2 but with some additional features:

- Masks bases below a minimum quality threshold before voting
- Produces consensus qualities capped at a configurable maximum (default Q90 for duplex)
- Requires minimum reads per strand (`min-reads` parameter)
- Handles insertions and deletions by alignment within the family

kam doesn't need to replicate fgbio exactly, but these design choices are worth noting as prior art.

## Duplex Crossing: Combining Forward and Reverse SSC

After computing SSC for each strand, the duplex consensus combines them. The key insight: if both strands independently call the same base, the probability that it's wrong is the product of both error probabilities — this is the ~10⁻⁷ figure for duplex sequencing.

### Agreement

```
SSC_fwd: G with error_prob 0.001
SSC_rev: G with error_prob 0.002

Duplex error_prob = SSC_fwd_error × SSC_rev_error = 0.001 × 0.002 = 0.000002
Duplex Phred = -10 × log10(0.000002) ≈ Q57
```

### Disagreement

```
SSC_fwd: G with error_prob 0.001
SSC_rev: A with error_prob 0.01
```

Options for handling this, in order of conservatism:

**Option A: Mask as N**
- Safest. Position is ambiguous, don't guess.
- Lost information — downstream k-mers spanning this position are destroyed.

**Option B: Pick higher-confidence SSC**
- Here: pick G (lower error probability from forward strand)
- Consensus quality reflects only the winning strand, NOT duplex-level confidence
- Mark as `SIMPLEX_RESOLVED` in per-base annotation so downstream knows this position doesn't have duplex support

**Option C: Report both**
- Emit the position with a flag indicating strand disagreement
- Let the variant caller handle it — a position with strand disagreement at the variant site is suspicious

**Recommendation:** Option B as default (pick higher-confidence), with the position flagged as non-duplex. Option A available via `--strict-duplex` flag. This preserves the most information while being honest about confidence.

## Handling Edge Cases

### Family Size = 1 (Singleton)

No consensus to compute — the single read IS the "consensus." Quality is just the raw read quality. Typically filtered out by minimum family size thresholds, but the parser should handle this gracefully.

### Family Size = 2 with Disagreement

```
Read 1: G at Q30
Read 2: A at Q30
```

50/50 split with equal quality. Both algorithms pick arbitrarily (or the base that comes first alphabetically). Consensus quality is low (~Q3). This is a genuinely ambiguous position.

**Recommendation:** Flag these positions explicitly. A family of 2 with a disagreement is not much better than a singleton. If `--min-family-size 3` is set, this family might be filtered anyway.

### Systematic Errors (Same Error in Multiple Reads)

If a sequence context causes the same error in multiple reads (e.g., a homopolymer run that consistently miscalls), majority vote will call the error as consensus. This is the fundamental limitation of single-strand consensus.

Duplex consensus catches this because the error is typically strand-specific — the forward strand has the systematic error, but the reverse strand doesn't. This is one of the strongest arguments for requiring duplex support for high-confidence calls.

### Insertions and Deletions Within a Family

Reads from the same molecule should have the same length, but indel errors can cause length variation. Options:
- **Simple:** require all reads in a family to have the same length. Reads with different lengths are flagged as potential alignment artefacts.
- **Complex:** perform multiple sequence alignment (MSA) within the family to handle indel errors. Significantly more computation.

**Recommendation for initial implementation:** require same-length reads within a family. Flag length mismatches in logs. MSA-based consensus as a later enhancement if needed.

## Recommendation for kam

### Phase 1: Quality-Weighted Majority Vote (Algorithm 2)
- Simple, fast, good enough for initial pipeline validation
- Easy to test with synthetic data (expected consensus is deterministic)
- Per-base consensus quality estimated from vote margin

### Phase 2: Bayesian Posterior (Algorithm 3)
- Replace Algorithm 2's quality estimation with proper Bayesian posteriors
- Produces calibrated confidence scores for downstream statistical calling
- Log-space computation for numerical stability

### Duplex Crossing: Option B (higher-confidence SSC wins at disagreements)
- Flag non-duplex positions in per-base annotation
- `--strict-duplex` flag for Option A (mask as N)
- Report strand disagreement rate in QC

### Configurable Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--consensus-algorithm` | `quality-weighted` | `quality-weighted` or `bayesian` |
| `--min-base-quality` | 10 | Mask bases below this quality before voting |
| `--max-consensus-quality` | 60 | Cap consensus Phred score (Illumina can't measure above ~Q40 anyway) |
| `--strict-duplex` | false | Mask disagreement positions as N instead of picking best strand |
| `--min-family-size` | 1 | Minimum reads to attempt consensus (1 = keep singletons) |
| `--min-duplex-reads` | 0 | Minimum reads per strand for duplex (0 = allow simplex) |

## References

Schmitt, M.W., Kennedy, S.R., Salk, J.J., Fox, E.J., Hiatt, J.B. and Loeb, L.A. (2012) 'Detection of ultra-rare mutations by next-generation sequencing', Proceedings of the National Academy of Sciences, 109(36), pp. 14508–14513. doi: 10.1073/pnas.1208715109. [Foundational duplex sequencing paper establishing the ~10⁻⁷ error rate for duplex consensus and the theoretical basis for strand-specific error cancellation.]

Fulcrum Genomics (2024) fgbio toolkit. GitHub repository. Available at: https://github.com/fulcrumgenomics/fgbio (Accessed: March 2026). [Source of the `CallDuplexConsensusReads` implementation details: quality masking, consensus quality capping, `min-reads` parameter design, and the duplex crossing algorithm that kam's approach is compared against.]

Kennedy, S.R., Schmitt, M.W., Fox, E.J., Kohrn, B.F., Salk, J.J., Ahn, E.H., Prindle, M.J., Kuong, K.J., Shen, J.-C., Risques, R.-A. and Loeb, L.A. (2014) 'Detecting ultralow-frequency mutations by Duplex Sequencing', Nature Protocols, 9(11), pp. 2586–2606. doi: 10.1038/nprot.2014.170. [Detailed protocol for duplex sequencing including the per-strand consensus → duplex consensus two-level approach and the mathematical basis for error probability multiplication across strands.]
