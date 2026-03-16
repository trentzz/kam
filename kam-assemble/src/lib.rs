//! Molecule assembly from raw FASTQ for Twist UMI duplex sequencing.
//!
//! `kam-assemble` replaces HUMID with full molecule-level evidence preservation.
//! It handles inline UMI extraction, canonical UMI pairing, Hamming-distance
//! clustering, and quality-weighted consensus calling.

/// Placeholder module — molecule assembly will be implemented in subsequent tasks.
pub mod assembly {}

pub mod assembler;
pub mod clustering;
pub mod consensus;
pub mod fingerprint;

/// FASTQ I/O for paired-end reads using needletail.
///
/// See [`io::read_fastq_pairs`] for the main entry point.
pub mod io;

/// Twist UMI FASTQ read pair parser.
///
/// See [`parser::parse_read_pair`] for the main entry point.
pub mod parser;
