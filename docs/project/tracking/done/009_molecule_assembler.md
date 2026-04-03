# Task 009: Molecule Assembly Orchestrator

## Location
`kam-assemble/src/assembler.rs`

## What to implement

The main orchestrator that ties together parser, clustering, fingerprinting, and consensus calling into a complete molecule assembly pipeline. Takes paired FASTQ reads, outputs assembled Molecules.

## Interface

```rust
use kam_core::molecule::{Molecule, ConsensusRead, CanonicalUmiPair, Strand, FamilyType};
use crate::parser::{ParserConfig, ParseResult, ParseStats, ParsedReadPair};
use crate::clustering::cluster_umi_pairs;
use crate::fingerprint::{compute_endpoint_fingerprint, fingerprints_compatible};
use crate::consensus::{ConsensusConfig, single_strand_consensus, duplex_consensus, DisagreementStrategy};

/// Full configuration for molecule assembly
#[derive(Debug, Clone)]
pub struct AssemblerConfig {
    pub parser: ParserConfig,
    pub consensus: ConsensusConfig,
    pub disagreement_strategy: DisagreementStrategy,
    pub max_hamming_distance: u32,       // for UMI clustering, default 1
    pub min_family_size: u8,             // minimum reads to form a family, default 1
    pub min_duplex_reads: u8,            // minimum reads per strand for duplex, default 1
}

/// Statistics from assembly run
#[derive(Debug, Clone, Default)]
pub struct AssemblyStats {
    pub parse_stats: ParseStats,
    pub n_molecules: u64,
    pub n_duplex: u64,
    pub n_simplex_fwd: u64,
    pub n_simplex_rev: u64,
    pub n_singletons: u64,
    pub n_families_below_min_size: u64,
    pub n_umi_collisions_detected: u64,
}

/// Assemble molecules from parsed read pairs.
/// Groups reads by canonical UMI (with Hamming distance clustering),
/// sub-groups by endpoint fingerprint (collision detection),
/// then calls consensus per family.
pub fn assemble_molecules(
    read_pairs: Vec<ParsedReadPair>,
    config: &AssemblerConfig,
) -> (Vec<Molecule>, AssemblyStats);
```

## Algorithm

1. Group `ParsedReadPair`s by `canonical_umi`
2. Within each UMI group, apply Hamming distance clustering (`cluster_umi_pairs`)
3. Within each cluster, sub-group by endpoint fingerprint compatibility (collision detection)
4. For each final group:
   a. Separate into forward and reverse strand reads
   b. If total reads < `min_family_size`, skip (count in stats)
   c. Call SSC for forward strand reads
   d. Call SSC for reverse strand reads
   e. If both strands have reads, call duplex consensus
   f. Classify family type (Duplex, SimplexFwd, SimplexRev, Singleton)
   g. Build `Molecule` struct

## Important

The `Molecule` struct in kam-core may need to be adjusted — check what fields it has and populate them from the assembly results. The `ConsensusRead` struct should be populated from the consensus calling output.

## Tests required

1. Single read pair → singleton molecule
2. Two read pairs with same UMI, same strand → simplex family with consensus
3. Two read pairs with same UMI, opposite strands → duplex molecule
4. Two read pairs with same UMI but different fingerprints → split into two molecules (collision detected)
5. UMI clustering: two pairs with UMIs at Hamming distance 1 → merged into one molecule
6. `min_family_size` filters small families
7. AssemblyStats counts are correct
8. Empty input returns empty vec

## Definition of done

- `cargo test -p kam-assemble` passes
- `cargo clippy -p kam-assemble -- -D warnings` passes
- Doc comments with examples
