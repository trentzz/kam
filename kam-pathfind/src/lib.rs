//! De Bruijn graph construction and path walking.
//!
//! `kam-pathfind` builds a de Bruijn graph from a k-mer index and walks paths
//! between anchor k-mers in targeted mode to identify variant sequences.

pub mod anchor;
pub mod graph;
pub mod score;
pub mod walk;
