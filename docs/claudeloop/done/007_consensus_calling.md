# Task 007: Quality-Weighted Consensus Calling

## Location
`kam-assemble/src/consensus.rs`

## What to implement

Single-strand consensus (SSC) calling: given a family of reads from one strand, produce a consensus sequence with per-base error probabilities using quality-weighted majority vote.

See `docs/research/consensus_calling_algorithms.md` for algorithm details.

## Interface

```rust
/// Per-base consensus result
#[derive(Debug, Clone)]
pub struct ConsensusBase {
    pub base: u8,              // winning base (A/C/G/T)
    pub error_prob: f32,       // estimated probability of error
    pub depth: u8,             // number of reads covering this position
    pub votes: [u32; 4],       // weighted votes for A, C, G, T (scaled to integer)
}

/// Configuration for consensus calling
#[derive(Debug, Clone)]
pub struct ConsensusConfig {
    pub min_base_quality: u8,      // mask bases below this Phred score (default: 10)
    pub max_consensus_quality: u8,  // cap output Phred (default: 60)
}

impl Default for ConsensusConfig {
    fn default() -> Self {
        Self {
            min_base_quality: 10,
            max_consensus_quality: 60,
        }
    }
}

/// Convert Phred quality score to error probability
pub fn phred_to_prob(phred: u8) -> f32;

/// Convert error probability to Phred quality score (capped at max)
pub fn prob_to_phred(prob: f32, max_phred: u8) -> u8;

/// Compute single-strand consensus from a set of sequences and qualities.
/// All sequences must be the same length.
/// Returns None if sequences is empty.
pub fn single_strand_consensus(
    sequences: &[&[u8]],
    qualities: &[&[u8]],       // Phred+33 encoded
    config: &ConsensusConfig,
) -> Option<Vec<ConsensusBase>>;
```

## Algorithm

At each position:
1. For each read, convert Phred quality to error probability: `p = 10^(-Q/10)`
2. Weight = `1 - p` for the called base
3. Sum weights per base (A, C, G, T)
4. Consensus base = argmax of summed weights
5. Consensus error prob = `1 - (winner_weight / total_weight)`
6. Skip bases below `min_base_quality` (don't include in voting)
7. Cap output quality at `max_consensus_quality`

## Tests required

1. Single read → consensus equals that read
2. Three reads agreeing → consensus matches, quality improves
3. Three reads, one disagreeing at low quality → consensus follows majority
4. Two reads disagreeing at equal quality → picks one, low confidence
5. `phred_to_prob(30)` ≈ 0.001
6. `phred_to_prob(20)` ≈ 0.01
7. `prob_to_phred(0.001, 60)` = 30
8. `prob_to_phred` caps at max_phred
9. Base below `min_base_quality` is excluded from voting
10. Empty input returns None
11. All bases at a position masked → ConsensusBase has error_prob = 1.0

## Definition of done

- `cargo test -p kam-assemble` passes
- `cargo clippy -p kam-assemble -- -D warnings` passes
- Doc comments with examples
