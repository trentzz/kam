# Correctness Review

## Algorithmic Correctness

### UMI Parsing and Canonical Pairing — Correct

The canonical UMI pairing logic is correct: `min(R1_UMI + R2_UMI, R2_UMI + R1_UMI)` lexicographically. The implementation uses byte-level comparison which matches lexicographic ordering for ASCII DNA bases. The `strand_of_r1()` function correctly identifies which strand R1 came from by checking which UMI is `umi_a`.

### Hamming Distance Clustering — Correct but Incomplete

The directional clustering algorithm is correct: higher-count seeds absorb lower-count neighbors within the Hamming distance threshold. No transitive chaining is allowed. This matches the UMI-tools directional method.

**Concern**: The algorithm compares combined 10-base UMI pairs. A Hamming distance of 1 on a 10-base pair is very conservative — it means only 1 base in 10 can differ. The research docs discuss this tradeoff but the default `max_hamming_distance = 1` may under-cluster in practice, especially at high PCR duplicate rates where sequencing errors in UMIs are common.

### Consensus Calling — Mostly Correct

The quality-weighted majority vote is correctly implemented. The weight `w = 1 - 10^(-Q/10)` is standard. The argmax selects the base with the highest total weight.

**Issue**: The error probability formula is questionable:

```rust
let raw_error = 1.0 - winner_weight / total_weight;
```

This says "error probability = fraction of weight that didn't go to the winner." For 3 reads of ACGT quality calling A, the total weight is `3 * 0.999 = 2.997`, the winner weight is `2.997`, and the error prob is `0.0`. That's correct for unanimous agreement.

For 2 reads calling A at Q30 and 1 read calling T at Q30: winner weight = `2 * 0.999 = 1.998`, total = `3 * 0.999 = 2.997`, error = `1 - 1.998/2.997 = 0.333`. This means 33% error probability, which seems too high for a 2:1 majority at Q30. The true Bayesian posterior would be much lower.

The formula treats quality weights as votes, not as probabilistic evidence. A proper Bayesian approach would multiply likelihoods across reads, not sum weights. This may matter for small families (2-3 reads) where the difference between vote-counting and Bayesian updating is significant.

### Duplex Consensus — Correct

Error probability multiplication for agreeing bases is correct: if fwd has error prob `p1` and rev has `p2`, the duplex error is `p1 * p2`. This is the standard independent-error model.

The `PickBest` disagreement strategy correctly selects the lower error probability. `MaskAsN` correctly sets error to 1.0.

**Issue**: `assert_eq!` is used to check fwd and rev SSC lengths match. This panics in library code. If templates have different lengths due to indels, this would crash. It should return an error.

### K-mer Encoding — Correct

The 2-bit encoding is standard (A=00, C=01, G=10, T=11). The sliding window iterator is O(1) per k-mer. The canonical form (`min(fwd, revcomp)`) is computed correctly. The reverse complement XOR trick (`encoded ^ mask`) correctly complements all bases, and the bit-reversal loop correctly reverses base order.

### De Bruijn Graph — Correct

Edge derivation uses `suffix(A, k-1) == prefix(B, k-1)` which is the standard de Bruijn graph construction. The prefix-grouping optimization avoids O(N^2) all-pairs comparison. Self-loops are handled correctly.

### BFS Path Walking — Correct with Caveat

The BFS with per-path visited sets correctly avoids cycles. The `max_path_length` and `max_paths` bounds prevent unbounded exploration.

**Caveat**: BFS finds *shortest paths first*, which is good for finding reference-length paths. But it may miss long insertion paths before hitting `max_paths`. If a region has many short paths (e.g., a repeat region), all 100 allowed paths might be short alternatives, and the biologically interesting long path might never be explored. DFS with depth limits or bidirectional BFS could be more appropriate.

### VAF Estimation — Correct

```rust
let point = k as f64 / m as f64;
let alpha = k as f64 + 1.0;
let beta_param = (m - k) as f64 + 1.0;
```

This is the conjugate Beta posterior with uniform prior `Beta(1, 1)` updated by `k` successes in `m` trials, giving `Beta(k+1, m-k+1)`. The 2.5% and 97.5% quantiles give the 95% credible interval. This is textbook correct.

**Issue**: The prior `Beta(1, 1)` is uninformative. For ctDNA variant detection at 0.1% VAF, a more informative prior (e.g., `Beta(0.5, 500)` reflecting the expectation that most sites are wildtype) would give tighter intervals and better calibration.

### Fisher's Exact Test — Correct Implementation

The hypergeometric log-PMF calculation is standard:

```rust
fn hypergeometric_log_pmf(a: i64, row1: i64, col1: i64, n: i64) -> f64 {
    log_binom(row1, a) + log_binom(n - row1, col1 - a) - log_binom(n, col1)
}
```

The two-tailed p-value is computed by summing all tables with probability ≤ the observed table's probability. This is the standard Clopper-Pearson approach.

**Issue**: `log_factorial(n)` uses a naive sum of logarithms: `(2..=n).map(|i| (i as f64).ln()).sum()`. For large n (say, 10,000 molecules), this is 10,000 log computations. Stirling's approximation or `lgamma` would be more efficient. For the typical panel use case (depth ~1000), this is fine, but it's O(n) per call.

### Confidence Scoring — Simplified

```rust
fn compute_confidence(k: u32, m: u32, background_error_rate: f64) -> f64 {
    let log_l_signal = log_binomial_likelihood(k, m, vaf);
    let log_l_background = log_binomial_likelihood(k, m, background_error_rate);
    let log_ratio = log_l_signal - log_l_background;
    exp_ratio / (1.0 + exp_ratio)
}
```

This is a binomial likelihood ratio test converted to a posterior probability with equal priors. It's mathematically correct but oversimplified:

1. The signal hypothesis uses the observed VAF as the true frequency, which overestimates the likelihood ratio (circular reasoning)
2. Equal priors (`P(signal) = P(noise) = 0.5`) are unrealistic — most sites are wildtype
3. No overdispersion modeling (beta-binomial would be more appropriate)
4. No PCR error vs sequencing error distinction

For a first implementation, this is acceptable. The docs correctly identify Phase 2 as "beta-binomial with per-site background."

### Variant Classification — Correct

```rust
match ref_seq.len().cmp(&alt_seq.len()) {
    Ordering::Equal => { /* count diffs: 0-1 = SNV, >1 = MNV */ }
    Ordering::Greater => VariantType::Deletion,
    Ordering::Less => VariantType::Insertion,
}
```

This is simple and correct. One edge: `diffs == 0` (identical sequences) is classified as SNV, which is wrong — it should be `None` or `Reference`. This would happen if the path walking returns a path identical to the reference.

### Endpoint Fingerprinting — Correct

The 2-bit encoding, 4-endpoint packing, and popcount-based compatibility check are all correct. The 4-bit threshold tolerates ~2 base errors across 32 endpoint bases, which is reasonable for Q20+ data.

**Edge case**: Templates shorter than 16bp have overlapping first/last endpoints. The code handles this correctly (the fingerprint has repeated information but doesn't panic). However, very short templates produce low-entropy fingerprints that may be falsely compatible with unrelated molecules.

## Biological Correctness Concerns

### 1. Forward/Reverse Strand Assignment is Arbitrary

The canonical UMI pairing assigns "Forward" to whichever read has the lexicographically smaller UMI. This has no biological meaning — it's purely a convention for consistent grouping. The code correctly uses this only for grouping and never interprets it as actual DNA strand orientation.

### 2. Template Orientation Not Reversed

When building duplex consensus, the code uses `template_r1` for both forward and reverse strand reads:

```rust
// kam-assemble/src/assembler.rs:443
let seqs: Vec<&[u8]> = reads.iter().map(|rp| rp.template_r1.as_slice()).collect();
```

For duplex sequencing, the reverse strand read should be reverse-complemented before alignment with the forward strand. If `template_r1` from a reverse-strand read is used as-is, the duplex consensus comparison is comparing sequences in opposite orientations, which would always disagree.

**This may be the most serious correctness issue in the codebase.** If templates are not reverse-complemented before duplex comparison, the duplex consensus will never agree, and duplex error correction (the key advantage of duplex sequencing) won't work.

However, it's possible that needletail or the read structure already handles this — if R1 and R2 are already in the same orientation for Twist chemistry. This needs careful verification against actual Twist UMI data.

### 3. Family Size Tracked as u8

```rust
pub family_size: (u8, u8),
```

`u8` limits family size to 255 reads per strand. In high-depth sequencing (e.g., deep ctDNA panels at 10,000x), families can exceed this. The parser also uses `u8`:

```rust
let n_fwd = fwd_reads.len() as u8;
```

This will silently truncate at 255. Use `u16` or `u32`.

### 4. No Read Name Tracking

The pipeline discards read names immediately. For debugging, QC, and regulatory compliance, the ability to trace a variant call back to specific read pairs is important. At minimum, a debug mode should log which reads contributed to each molecule.
