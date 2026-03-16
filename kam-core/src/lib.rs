//! Shared types and traits for the kam alignment-free variant detection pipeline.
//!
//! `kam-core` defines the fundamental data model that all other crates build on.
//! The [`molecule::Molecule`] is the atomic unit of evidence throughout the
//! pipeline.
//!
//! ## Crate layout
//!
//! | Module | Contents |
//! |---|---|
//! | [`molecule`] | [`molecule::CanonicalUmiPair`], [`molecule::Strand`], [`molecule::FamilyType`], [`molecule::ConsensusRead`] |
//! | [`kmer`] | [`kmer::MoleculeEvidence`], [`kmer::KmerIndex`] trait |
//! | [`chemistry`] | [`chemistry::ReadStructure`] with Twist UMI duplex preset |
//! | [`error`] | [`error::KamError`] pipeline error enum |
//! | [`qc`] | [`qc::AssemblyQc`], [`qc::IndexQc`], [`qc::PathfindQc`], [`qc::CallQc`], [`qc::write_qc`] |
//! | [`serialize`] | [`serialize::FileType`], [`serialize::FileHeader`], [`serialize::write_bincode`], [`serialize::read_header`], [`serialize::read_bincode`] |

pub mod chemistry;
pub mod error;
pub mod kmer;
pub mod molecule;
pub mod qc;
pub mod serialize;
