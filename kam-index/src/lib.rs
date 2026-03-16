//! K-mer indexing with molecule provenance.
//!
//! `kam-index` replaces Jellyfish for this use case, storing `MoleculeEvidence`
//! per k-mer rather than raw integer counts. Supports allowlist filtering
//! to bound memory usage to panel size rather than sequencing depth.

pub mod allowlist;
pub mod encode;
pub mod extract;
pub mod hash_index;

pub use hash_index::HashKmerIndex;
