# Task 008: Duplex Consensus Crossing

## Location
`kam-assemble/src/consensus.rs` (extend existing module)

## What to implement

Duplex consensus: combine forward-strand SSC and reverse-strand SSC into a duplex consensus. Positions where both strands agree get dramatically reduced error probability. Positions where they disagree are handled according to configuration.

## Interface

```rust
/// How to handle positions where forward and reverse SSC disagree
#[derive(Debug, Clone, Copy, Default)]
pub enum DisagreementStrategy {
    /// Pick the base from whichever SSC had lower error probability.
    /// Mark as non-duplex in output.
    #[default]
    PickBest,
    /// Mask as N with error_prob = 1.0
    MaskAsN,
}

/// Result of duplex consensus at a single position
#[derive(Debug, Clone)]
pub struct DuplexBase {
    pub base: u8,
    pub error_prob: f32,
    pub is_duplex: bool,        // true if both strands agreed
    pub fwd_depth: u8,
    pub rev_depth: u8,
}

/// Combine forward and reverse single-strand consensus into duplex consensus.
/// Both SSCs must be the same length.
/// Returns None if either input is empty.
pub fn duplex_consensus(
    fwd_ssc: &[ConsensusBase],
    rev_ssc: &[ConsensusBase],
    strategy: DisagreementStrategy,
) -> Option<Vec<DuplexBase>>;
```

## Algorithm

At each position:
1. If fwd and rev agree on the base:
   - `error_prob = fwd.error_prob * rev.error_prob` (independent errors multiply)
   - `is_duplex = true`
2. If fwd and rev disagree:
   - `PickBest`: use the base with lower error_prob, `is_duplex = false`, error_prob from winning strand only
   - `MaskAsN`: base = N, error_prob = 1.0, `is_duplex = false`

## Tests required

1. Both strands agree everywhere → all positions duplex, error probs multiplied
2. One position disagrees, PickBest → non-duplex at that position, uses better strand
3. One position disagrees, MaskAsN → N at that position, error_prob = 1.0
4. Error probability multiplication: 0.001 * 0.001 = 0.000001
5. Empty input returns None
6. Mismatched lengths panics or returns error
7. fwd_depth and rev_depth propagated correctly from SSCs

## Definition of done

- `cargo test -p kam-assemble` passes
- `cargo clippy -p kam-assemble -- -D warnings` passes
- Doc comments with examples
