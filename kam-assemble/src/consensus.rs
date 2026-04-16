//! Quality-weighted consensus calling and duplex consensus crossing.
//!
//! This module implements two layers of consensus:
//!
//! 1. **Single-strand consensus (SSC)** — [`single_strand_consensus`] collapses a
//!    family of reads from one strand into a per-base consensus using
//!    quality-weighted majority vote.
//!
//! 2. **Duplex consensus** — [`duplex_consensus`] combines a forward-strand SSC and
//!    a reverse-strand SSC into a final duplex consensus. Positions where both strands
//!    agree receive dramatically reduced error probabilities (independent error
//!    multiplication). Disagreements are resolved according to [`DisagreementStrategy`].

// ── Single-strand consensus ──────────────────────────────────────────────────

/// Per-base result from single-strand consensus calling.
///
/// # Example
/// ```
/// use kam_assemble::consensus::{single_strand_consensus, ConsensusConfig};
///
/// let seqs  = &[b"ACGT".as_ref()];
/// let quals = &[b"IIII".as_ref()]; // Phred+33, Q=40
/// let result = single_strand_consensus(seqs, quals, &ConsensusConfig::default()).unwrap();
/// assert_eq!(result[0].base, b'A');
/// ```
#[derive(Debug, Clone)]
pub struct ConsensusBase {
    /// Winning base (one of `A`, `C`, `G`, `T`, or `N` when all positions are masked).
    pub base: u8,
    /// Estimated probability of error at this position (0.0–1.0).
    pub error_prob: f32,
    /// Number of reads covering this position (including masked bases).
    pub depth: u8,
    /// Quality-weighted votes for `[A, C, G, T]` (scaled to integer, thousandths).
    pub votes: [u32; 4],
}

/// Configuration for consensus calling.
///
/// # Example
/// ```
/// use kam_assemble::consensus::ConsensusConfig;
///
/// let cfg = ConsensusConfig::default();
/// assert_eq!(cfg.min_base_quality, 10);
/// assert_eq!(cfg.max_consensus_quality, 60);
/// ```
#[derive(Debug, Clone)]
pub struct ConsensusConfig {
    /// Bases with Phred quality below this threshold are excluded from voting (default: 10).
    pub min_base_quality: u8,
    /// Output consensus Phred quality is capped at this value (default: 60).
    pub max_consensus_quality: u8,
}

impl Default for ConsensusConfig {
    fn default() -> Self {
        Self {
            min_base_quality: 10,
            max_consensus_quality: 60,
        }
    }
}

/// Convert a Phred quality score to an error probability.
///
/// `p = 10^(-Q/10)`
///
/// # Example
/// ```
/// use kam_assemble::consensus::phred_to_prob;
///
/// let p = phred_to_prob(30);
/// assert!((p - 0.001).abs() < 1e-6);
/// ```
pub fn phred_to_prob(phred: u8) -> f32 {
    10_f32.powf(-(phred as f32) / 10.0)
}

/// Convert an error probability to a Phred quality score, capped at `max_phred`.
///
/// `Q = -10 * log10(p)`, then clamped to `[0, max_phred]`.
///
/// # Example
/// ```
/// use kam_assemble::consensus::prob_to_phred;
///
/// assert_eq!(prob_to_phred(0.001, 60), 30);
/// assert_eq!(prob_to_phred(1e-10, 40), 40); // capped
/// ```
pub fn prob_to_phred(prob: f32, max_phred: u8) -> u8 {
    if prob <= 0.0 {
        return max_phred;
    }
    let q = (-10.0 * prob.log10()).round() as u8;
    q.min(max_phred)
}

/// Compute single-strand consensus from a set of sequences and qualities.
///
/// All sequences (and quality arrays) must be the same length. Bases with
/// Phred quality below `config.min_base_quality` are excluded from voting.
///
/// Returns `None` if `sequences` is empty.
///
/// # Algorithm
///
/// At each position:
/// 1. Convert each base's Phred quality to error probability `p = 10^(-Q/10)`.
/// 2. Weight for the called base = `1 - p`.
/// 3. Sum weights per nucleotide (A/C/G/T).
/// 4. Consensus base = argmax of summed weights.
/// 5. Consensus error prob = `1 - winner_weight / total_weight`.
///    If all bases at a position are masked the error prob is 1.0 and base is `N`.
/// 6. Output quality is capped at `config.max_consensus_quality`.
///
/// # Example
/// ```
/// use kam_assemble::consensus::{single_strand_consensus, ConsensusConfig};
///
/// let seqs  = &[b"ACGT".as_ref(), b"ACGT".as_ref()];
/// let quals = &[b"IIII".as_ref(), b"IIII".as_ref()]; // Q=40
/// let result = single_strand_consensus(seqs, quals, &ConsensusConfig::default()).unwrap();
/// assert_eq!(result[0].base, b'A');
/// assert!(result[0].error_prob < 0.001);
/// ```
pub fn single_strand_consensus(
    sequences: &[&[u8]],
    qualities: &[&[u8]],
    config: &ConsensusConfig,
) -> Option<Vec<ConsensusBase>> {
    if sequences.is_empty() {
        return None;
    }

    let seq_len = sequences[0].len();
    let mut result = Vec::with_capacity(seq_len);

    for pos in 0..seq_len {
        // Accumulated weighted votes for A=0, C=1, G=2, T=3.
        // Stored as f32 internally; exported as scaled integer.
        let mut weight_sum = [0_f32; 4];
        let mut total_weight = 0_f32;
        let depth = sequences.len().min(u8::MAX as usize) as u8;

        for (seq, qual) in sequences.iter().zip(qualities.iter()) {
            let base = seq[pos];
            let q_raw = qual[pos].saturating_sub(33); // Phred+33 → Phred
            if q_raw < config.min_base_quality {
                continue;
            }
            let p = phred_to_prob(q_raw);
            let w = 1.0 - p;
            let idx = base_to_idx(base);
            if let Some(i) = idx {
                weight_sum[i] += w;
                total_weight += w;
            }
        }

        let (consensus_base, error_prob) = if total_weight <= 0.0 {
            // All bases masked.
            (b'N', 1.0_f32)
        } else {
            // argmax
            let (winner_idx, winner_weight) = weight_sum
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(i, &w)| (i, w))
                .expect("weight_sum is non-empty when total_weight > 0");
            let raw_error = 1.0 - winner_weight / total_weight;
            // Floor at 0 to avoid -0.0 artefacts from floating point.
            let error_prob = raw_error.max(0.0);
            (idx_to_base(winner_idx), error_prob)
        };

        // Cap output quality.
        let capped_error = if error_prob < 1.0 {
            let phred = prob_to_phred(error_prob, config.max_consensus_quality);
            phred_to_prob(phred)
        } else {
            1.0
        };

        // Scale weights to integer (thousandths of a vote unit).
        let votes = [
            (weight_sum[0] * 1000.0).round() as u32,
            (weight_sum[1] * 1000.0).round() as u32,
            (weight_sum[2] * 1000.0).round() as u32,
            (weight_sum[3] * 1000.0).round() as u32,
        ];

        result.push(ConsensusBase {
            base: consensus_base,
            error_prob: capped_error,
            depth,
            votes,
        });
    }

    Some(result)
}

// ── Helpers ──────────────────────────────────────────────────────────────────

fn base_to_idx(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn idx_to_base(idx: usize) -> u8 {
    match idx {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'N',
    }
}

// ── Duplex consensus ─────────────────────────────────────────────────────────

/// How to handle positions where forward and reverse SSC disagree.
///
/// # Example
/// ```
/// use kam_assemble::consensus::DisagreementStrategy;
///
/// let strategy = DisagreementStrategy::default(); // PickBest
/// let _ = strategy; // used by duplex_consensus
/// ```
#[derive(Debug, Clone, Copy, Default)]
pub enum DisagreementStrategy {
    /// Pick the base from whichever SSC had lower error probability.
    /// The position is marked as non-duplex.
    #[default]
    PickBest,
    /// Mask the position as `N` with `error_prob = 1.0`.
    /// The position is marked as non-duplex.
    MaskAsN,
}

/// Per-base result of duplex consensus calling.
///
/// # Example
/// ```
/// use kam_assemble::consensus::{
///     DuplexBase, ConsensusBase, DisagreementStrategy, duplex_consensus,
/// };
/// ```
#[derive(Debug, Clone)]
pub struct DuplexBase {
    /// Called base.
    pub base: u8,
    /// Estimated error probability (0.0–1.0).
    pub error_prob: f32,
    /// `true` when both strands agreed on this position.
    pub is_duplex: bool,
    /// Depth from the forward-strand SSC.
    pub fwd_depth: u8,
    /// Depth from the reverse-strand SSC.
    pub rev_depth: u8,
}

/// Combine forward and reverse single-strand consensus into a duplex consensus.
///
/// Both SSCs must be the same length. Returns `None` if either input slice is
/// empty or if the two SSCs have different lengths.
///
/// # Algorithm
///
/// At each position:
/// - If both strands agree: `error_prob = fwd.error_prob * rev.error_prob`, `is_duplex = true`.
/// - If they disagree:
///   - [`DisagreementStrategy::PickBest`]: use the base with lower `error_prob`, `is_duplex = false`.
///   - [`DisagreementStrategy::MaskAsN`]: base = `N`, `error_prob = 1.0`, `is_duplex = false`.
///
/// # Example
/// ```
/// use kam_assemble::consensus::{
///     single_strand_consensus, duplex_consensus, ConsensusConfig, DisagreementStrategy,
/// };
///
/// let seqs  = &[b"ACGT".as_ref()];
/// let quals = &[b"IIII".as_ref()];
/// let cfg = ConsensusConfig::default();
/// let fwd = single_strand_consensus(seqs, quals, &cfg).unwrap();
/// let rev = single_strand_consensus(seqs, quals, &cfg).unwrap();
/// let duplex = duplex_consensus(&fwd, &rev, DisagreementStrategy::PickBest).unwrap();
/// assert!(duplex[0].is_duplex);
/// ```
pub fn duplex_consensus(
    fwd_ssc: &[ConsensusBase],
    rev_ssc: &[ConsensusBase],
    strategy: DisagreementStrategy,
) -> Option<Vec<DuplexBase>> {
    if fwd_ssc.is_empty() || rev_ssc.is_empty() {
        return None;
    }

    // A length mismatch is unexpected (both SSCs should derive from the same
    // template), but must not panic in library code. Return None with a
    // diagnostic message so the caller can log and skip the molecule.
    if fwd_ssc.len() != rev_ssc.len() {
        log::warn!(
            "duplex_consensus: length mismatch (fwd={}, rev={}) — skipping molecule",
            fwd_ssc.len(),
            rev_ssc.len()
        );
        return None;
    }

    let result = fwd_ssc
        .iter()
        .zip(rev_ssc.iter())
        .map(|(fwd, rev)| {
            let fwd_depth = fwd.depth;
            let rev_depth = rev.depth;

            if fwd.base == rev.base {
                // Both strands agree — multiply independent error probabilities.
                DuplexBase {
                    base: fwd.base,
                    error_prob: fwd.error_prob * rev.error_prob,
                    is_duplex: true,
                    fwd_depth,
                    rev_depth,
                }
            } else {
                // Disagreement — apply chosen strategy.
                match strategy {
                    DisagreementStrategy::PickBest => {
                        let (base, error_prob) = if fwd.error_prob <= rev.error_prob {
                            (fwd.base, fwd.error_prob)
                        } else {
                            (rev.base, rev.error_prob)
                        };
                        DuplexBase {
                            base,
                            error_prob,
                            is_duplex: false,
                            fwd_depth,
                            rev_depth,
                        }
                    }
                    DisagreementStrategy::MaskAsN => DuplexBase {
                        base: b'N',
                        error_prob: 1.0,
                        is_duplex: false,
                        fwd_depth,
                        rev_depth,
                    },
                }
            }
        })
        .collect();

    Some(result)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Task 007: phred/prob conversions ────────────────────────────────────

    #[test]
    fn phred30_to_prob() {
        let p = phred_to_prob(30);
        assert!((p - 0.001).abs() < 1e-6, "expected ~0.001, got {p}");
    }

    #[test]
    fn phred20_to_prob() {
        let p = phred_to_prob(20);
        assert!((p - 0.01).abs() < 1e-5, "expected ~0.01, got {p}");
    }

    #[test]
    fn prob_to_phred_30() {
        assert_eq!(prob_to_phred(0.001, 60), 30);
    }

    #[test]
    fn prob_to_phred_capped() {
        // probability so small that the raw Phred would be >40
        assert_eq!(prob_to_phred(1e-10_f32, 40), 40);
    }

    // ── Task 007: single_strand_consensus ───────────────────────────────────

    #[test]
    fn single_read_identity() {
        let cfg = ConsensusConfig::default();
        let seq = b"ACGT";
        // Phred+33: '!' = 0, 'I' = 40
        let qual = b"IIII";
        let result = single_strand_consensus(&[seq.as_ref()], &[qual.as_ref()], &cfg).unwrap();
        assert_eq!(result.len(), 4);
        assert_eq!(result[0].base, b'A');
        assert_eq!(result[1].base, b'C');
        assert_eq!(result[2].base, b'G');
        assert_eq!(result[3].base, b'T');
    }

    #[test]
    fn three_agreeing_reads_improve_quality() {
        let cfg = ConsensusConfig::default();
        let seq = b"AAAA";
        // Q=30 → p=0.001 each
        let qual = b"????"; // '?' = Phred+33 → Q=30
        let seqs = &[seq.as_ref(); 3];
        let quals = &[qual.as_ref(); 3];
        let result = single_strand_consensus(seqs, quals, &cfg).unwrap();
        // Three agreeing reads: winner_weight = 3*(1-0.001)=2.997, total=2.997
        // error_prob ~ 0 → quality should be at or near max.
        for base_res in &result {
            assert_eq!(base_res.base, b'A');
            assert!(
                base_res.error_prob < 0.001,
                "error_prob={}",
                base_res.error_prob
            );
        }
    }

    #[test]
    fn majority_wins_over_low_quality_disagreement() {
        let cfg = ConsensusConfig::default();
        // Two reads say 'A' at high quality; one says 'T' at low quality.
        let seqs = &[b"A".as_ref(), b"A".as_ref(), b"T".as_ref()];
        // '?' = Q30, '+' = Phred+33 → Q=10 (just at threshold), '#' = Q=2 (below threshold)
        let quals = &[b"?".as_ref(), b"?".as_ref(), b"+".as_ref()];
        let result = single_strand_consensus(seqs, quals, &cfg).unwrap();
        assert_eq!(result[0].base, b'A');
    }

    #[test]
    fn two_equal_quality_reads_disagree_low_confidence() {
        let cfg = ConsensusConfig::default();
        // One 'A' and one 'T' at equal quality — picks one but confidence is low.
        let seqs = &[b"A".as_ref(), b"T".as_ref()];
        let quals = &[b"?".as_ref(), b"?".as_ref()]; // same Q=30
        let result = single_strand_consensus(seqs, quals, &cfg).unwrap();
        // Winner gets 50 % of weight → error_prob = 0.5
        assert!(
            result[0].error_prob >= 0.4 && result[0].error_prob <= 0.6,
            "expected ~0.5, got {}",
            result[0].error_prob
        );
    }

    #[test]
    fn base_below_min_quality_excluded() {
        let cfg = ConsensusConfig {
            min_base_quality: 10,
            max_consensus_quality: 60,
        };
        // One 'A' at Q=30; one 'T' at Q=2 (below threshold) — T vote must be excluded.
        let seqs = &[b"A".as_ref(), b"T".as_ref()];
        // '#' = Q=2 (Phred+33: '#' is ASCII 35, 35-33=2)
        let quals = &[b"?".as_ref(), b"#".as_ref()];
        let result = single_strand_consensus(seqs, quals, &cfg).unwrap();
        assert_eq!(result[0].base, b'A');
    }

    #[test]
    fn empty_input_returns_none() {
        let cfg = ConsensusConfig::default();
        let result = single_strand_consensus(&[], &[], &cfg);
        assert!(result.is_none());
    }

    #[test]
    fn all_bases_masked_gives_n_and_error_one() {
        let cfg = ConsensusConfig {
            min_base_quality: 20,
            max_consensus_quality: 60,
        };
        // Both reads have quality Q=2 (below threshold of 20).
        let seqs = &[b"A".as_ref(), b"A".as_ref()];
        let quals = &[b"#".as_ref(), b"#".as_ref()]; // '#' = Q=2
        let result = single_strand_consensus(seqs, quals, &cfg).unwrap();
        assert_eq!(result[0].base, b'N');
        assert_eq!(result[0].error_prob, 1.0);
    }

    // ── Task 008: duplex_consensus ───────────────────────────────────────────

    #[test]
    fn duplex_all_agree_is_duplex() {
        let cfg = ConsensusConfig::default();
        let seqs = &[b"ACGT".as_ref()];
        let quals = &[b"IIII".as_ref()]; // Q=40
        let fwd = single_strand_consensus(seqs, quals, &cfg).unwrap();
        let rev = single_strand_consensus(seqs, quals, &cfg).unwrap();
        let duplex = duplex_consensus(&fwd, &rev, DisagreementStrategy::PickBest).unwrap();
        assert_eq!(duplex.len(), 4);
        for d in &duplex {
            assert!(d.is_duplex, "expected all duplex");
        }
    }

    #[test]
    fn duplex_error_prob_multiplied() {
        // Build two ConsensusBase vectors with known error_prob = 0.001 each.
        let make_base = |ep: f32| ConsensusBase {
            base: b'A',
            error_prob: ep,
            depth: 3,
            votes: [3000, 0, 0, 0],
        };
        let fwd = vec![make_base(0.001)];
        let rev = vec![make_base(0.001)];
        let duplex = duplex_consensus(&fwd, &rev, DisagreementStrategy::PickBest).unwrap();
        let expected = 0.001_f32 * 0.001_f32;
        assert!(
            (duplex[0].error_prob - expected).abs() < 1e-12,
            "expected {expected}, got {}",
            duplex[0].error_prob
        );
    }

    #[test]
    fn duplex_disagreement_pick_best() {
        let fwd_base = ConsensusBase {
            base: b'A',
            error_prob: 0.01,
            depth: 2,
            votes: [2000, 0, 0, 0],
        };
        let rev_base = ConsensusBase {
            base: b'T',
            error_prob: 0.05,
            depth: 2,
            votes: [0, 0, 0, 2000],
        };
        let fwd = vec![fwd_base];
        let rev = vec![rev_base];
        let duplex = duplex_consensus(&fwd, &rev, DisagreementStrategy::PickBest).unwrap();
        assert_eq!(duplex[0].base, b'A'); // fwd has lower error
        assert_eq!(duplex[0].error_prob, 0.01);
        assert!(!duplex[0].is_duplex);
    }

    #[test]
    fn duplex_disagreement_mask_as_n() {
        let fwd_base = ConsensusBase {
            base: b'A',
            error_prob: 0.01,
            depth: 2,
            votes: [2000, 0, 0, 0],
        };
        let rev_base = ConsensusBase {
            base: b'T',
            error_prob: 0.05,
            depth: 2,
            votes: [0, 0, 0, 2000],
        };
        let fwd = vec![fwd_base];
        let rev = vec![rev_base];
        let duplex = duplex_consensus(&fwd, &rev, DisagreementStrategy::MaskAsN).unwrap();
        assert_eq!(duplex[0].base, b'N');
        assert_eq!(duplex[0].error_prob, 1.0);
        assert!(!duplex[0].is_duplex);
    }

    #[test]
    fn duplex_empty_returns_none() {
        let fwd: Vec<ConsensusBase> = vec![];
        let rev: Vec<ConsensusBase> = vec![];
        assert!(duplex_consensus(&fwd, &rev, DisagreementStrategy::PickBest).is_none());
    }

    // BUG-005: mismatched-length SSC inputs must not panic; the function must
    // return None and emit a diagnostic message instead.
    #[test]
    fn duplex_mismatched_lengths_returns_none() {
        let base = ConsensusBase {
            base: b'A',
            error_prob: 0.001,
            depth: 1,
            votes: [1000, 0, 0, 0],
        };
        let fwd = vec![base.clone(), base.clone()];
        let rev = vec![base];
        let result = duplex_consensus(&fwd, &rev, DisagreementStrategy::PickBest);
        assert!(
            result.is_none(),
            "mismatched lengths should return None, not panic"
        );
    }

    #[test]
    fn duplex_depths_propagated() {
        let fwd_base = ConsensusBase {
            base: b'A',
            error_prob: 0.001,
            depth: 5,
            votes: [5000, 0, 0, 0],
        };
        let rev_base = ConsensusBase {
            base: b'A',
            error_prob: 0.001,
            depth: 3,
            votes: [3000, 0, 0, 0],
        };
        let duplex =
            duplex_consensus(&[fwd_base], &[rev_base], DisagreementStrategy::PickBest).unwrap();
        assert_eq!(duplex[0].fwd_depth, 5);
        assert_eq!(duplex[0].rev_depth, 3);
    }

    // ── New edge-case tests ──────────────────────────────────────────────

    // Majority vote: 3 reads say A, 1 says C at equal quality. A wins.
    #[test]
    fn majority_vote_three_a_one_c() {
        let cfg = ConsensusConfig::default();
        let seqs = &[b"A".as_ref(), b"A".as_ref(), b"A".as_ref(), b"C".as_ref()];
        let quals = &[b"?".as_ref(), b"?".as_ref(), b"?".as_ref(), b"?".as_ref()]; // all Q=30
        let result = single_strand_consensus(seqs, quals, &cfg).unwrap();
        assert_eq!(result[0].base, b'A', "majority of 3 A vs 1 C should yield A");
        assert!(
            result[0].error_prob < 0.3,
            "error prob should be < 0.3 with 3:1 majority, got {}",
            result[0].error_prob
        );
    }

    // Tie-breaking is deterministic: 2 A vs 2 T at equal quality.
    // Run it 10 times to verify determinism.
    #[test]
    fn tie_breaking_deterministic() {
        let cfg = ConsensusConfig::default();
        let seqs = &[b"A".as_ref(), b"A".as_ref(), b"T".as_ref(), b"T".as_ref()];
        let quals = &[b"?".as_ref(), b"?".as_ref(), b"?".as_ref(), b"?".as_ref()];

        let results: Vec<u8> = (0..10)
            .map(|_| {
                single_strand_consensus(seqs, quals, &cfg)
                    .unwrap()[0]
                    .base
            })
            .collect();

        assert!(
            results.windows(2).all(|w| w[0] == w[1]),
            "tie-breaking must be deterministic: got {:?}",
            results
        );
    }

    // All bases are N at a position: result is N with error_prob = 1.0.
    #[test]
    fn all_n_bases_gives_n_and_error_one() {
        let cfg = ConsensusConfig::default();
        let seqs = &[b"N".as_ref(), b"N".as_ref(), b"N".as_ref()];
        let quals = &[b"I".as_ref(), b"I".as_ref(), b"I".as_ref()]; // Q=40
        let result = single_strand_consensus(seqs, quals, &cfg).unwrap();
        assert_eq!(result[0].base, b'N');
        assert_eq!(result[0].error_prob, 1.0, "all N bases should give error_prob = 1.0");
    }

    // Quality-weighted: a high-quality C beats a low-quality A.
    // A at Q=3 (w=0.5), C at Q=40 (w=0.9999). C wins.
    #[test]
    fn quality_weighted_high_quality_c_beats_low_quality_a() {
        let cfg = ConsensusConfig {
            min_base_quality: 1,
            max_consensus_quality: 60,
        };
        let seqs = &[b"A".as_ref(), b"C".as_ref()];
        // Q=3 for A (Phred+33 = ASCII 36 = '$'), Q=40 for C (Phred+33 = 'I')
        let quals = &[b"$".as_ref(), b"I".as_ref()];
        let result = single_strand_consensus(seqs, quals, &cfg).unwrap();
        assert_eq!(
            result[0].base, b'C',
            "high-quality C (Q=40) should beat low-quality A (Q=3)"
        );
    }

    // Depth field tracks total reads, including masked ones.
    #[test]
    fn depth_includes_masked_reads() {
        let cfg = ConsensusConfig {
            min_base_quality: 30,
            max_consensus_quality: 60,
        };
        let seqs = &[b"A".as_ref(), b"C".as_ref(), b"G".as_ref()];
        let quals = &[b"?".as_ref(), b"#".as_ref(), b"#".as_ref()]; // Q=30, Q=2, Q=2
        let result = single_strand_consensus(seqs, quals, &cfg).unwrap();
        assert_eq!(result[0].depth, 3, "depth should count all reads, including masked");
        assert_eq!(result[0].base, b'A', "only the A read passes quality threshold");
    }

    // Phred 0 conversion: probability should be 1.0.
    #[test]
    fn phred_zero_gives_prob_one() {
        let p = phred_to_prob(0);
        assert!((p - 1.0).abs() < 1e-6, "Q=0 should give p=1.0, got {p}");
    }

    // prob_to_phred with prob=0.0: returns max_phred.
    #[test]
    fn prob_zero_gives_max_phred() {
        assert_eq!(prob_to_phred(0.0, 60), 60);
        assert_eq!(prob_to_phred(0.0, 40), 40);
    }

    // prob_to_phred with prob=1.0: returns 0.
    #[test]
    fn prob_one_gives_phred_zero() {
        assert_eq!(prob_to_phred(1.0, 60), 0);
    }

    // Duplex agreement: error probability is product with different values.
    #[test]
    fn duplex_agreement_different_error_probs() {
        let fwd = vec![ConsensusBase {
            base: b'G',
            error_prob: 0.01,
            depth: 5,
            votes: [0, 0, 5000, 0],
        }];
        let rev = vec![ConsensusBase {
            base: b'G',
            error_prob: 0.001,
            depth: 3,
            votes: [0, 0, 3000, 0],
        }];
        let duplex = duplex_consensus(&fwd, &rev, DisagreementStrategy::PickBest).unwrap();
        let expected = 0.01_f32 * 0.001_f32;
        assert!(
            (duplex[0].error_prob - expected).abs() < 1e-10,
            "expected {expected}, got {}",
            duplex[0].error_prob
        );
        assert!(duplex[0].is_duplex);
    }

    // Duplex disagreement PickBest: when error_probs are equal, forward strand
    // is chosen (fwd.error_prob <= rev.error_prob).
    #[test]
    fn duplex_disagreement_pick_best_equal_error_picks_fwd() {
        let fwd = vec![ConsensusBase {
            base: b'A',
            error_prob: 0.05,
            depth: 2,
            votes: [2000, 0, 0, 0],
        }];
        let rev = vec![ConsensusBase {
            base: b'T',
            error_prob: 0.05,
            depth: 2,
            votes: [0, 0, 0, 2000],
        }];
        let duplex = duplex_consensus(&fwd, &rev, DisagreementStrategy::PickBest).unwrap();
        assert_eq!(
            duplex[0].base, b'A',
            "when error_probs are equal, PickBest should choose forward strand"
        );
        assert!(!duplex[0].is_duplex);
    }

    // Consensus with a single read: result is that read's sequence.
    #[test]
    fn single_read_consensus_matches_input() {
        let cfg = ConsensusConfig::default();
        let seq = b"GATTACA";
        let qual = b"IIIIIII"; // Q=40
        let result = single_strand_consensus(&[seq.as_ref()], &[qual.as_ref()], &cfg).unwrap();
        let consensus_seq: Vec<u8> = result.iter().map(|b| b.base).collect();
        assert_eq!(consensus_seq, b"GATTACA".to_vec());
    }

    // Multi-position consensus: each position resolved independently.
    #[test]
    fn multi_position_consensus_independent() {
        let cfg = ConsensusConfig::default();
        let seqs = &[b"AC".as_ref(), b"AC".as_ref(), b"TG".as_ref()];
        let quals = &[b"??".as_ref(), b"??".as_ref(), b"??".as_ref()];
        let result = single_strand_consensus(seqs, quals, &cfg).unwrap();
        assert_eq!(result[0].base, b'A');
        assert_eq!(result[1].base, b'C');
    }
}
