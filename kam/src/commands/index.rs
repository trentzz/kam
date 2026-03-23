//! Implementation of the `kam index` subcommand.
//!
//! Reads molecules from a bincode file, reads target sequences from FASTA,
//! builds an allowlist-filtered k-mer index, and writes the index to bincode.

use serde::{Deserialize, Serialize};

use kam_core::kmer::{KmerIndex, MoleculeEvidence};
use kam_core::molecule::{FamilyType, Molecule};
use kam_core::qc::{write_qc, IndexQc};
use kam_core::serialize::{read_bincode, write_bincode, FileType};
use kam_index::allowlist::build_allowlist;
use kam_index::extract::{extract_all, ConsensusReadInfo};
use kam_index::HashKmerIndex;

use crate::cli::IndexArgs;

/// A single entry in the serialized k-mer index file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerEntry {
    /// 2-bit encoded canonical k-mer.
    pub kmer: u64,
    /// Total molecules carrying this k-mer.
    pub n_molecules: u32,
    /// Duplex molecule count.
    pub n_duplex: u32,
    /// Simplex forward count.
    pub n_simplex_fwd: u32,
    /// Simplex reverse count.
    pub n_simplex_rev: u32,
    /// Minimum per-base error probability.
    pub min_base_error_prob: f32,
    /// Mean per-base error probability.
    pub mean_base_error_prob: f32,
}

/// Run the `index` subcommand end-to-end.
///
/// # Errors
///
/// Returns an error if file I/O or serialization fails.
pub fn run_index(args: IndexArgs) -> Result<(), Box<dyn std::error::Error>> {
    let k = args.kmer_size as usize;

    // ── 1. Read molecules from bincode ────────────────────────────────────────
    let (_header, molecules): (_, Vec<Molecule>) = read_bincode(&args.input)?;

    // ── 2. Read target sequences from FASTA ───────────────────────────────────
    let targets = read_fasta(&args.targets)?;

    // ── 3. Build allowlist from targets ───────────────────────────────────────
    let target_slices: Vec<&[u8]> = targets.iter().map(|(_id, seq)| seq.as_slice()).collect();
    let mut allowlist = build_allowlist(&target_slices, k);

    // Augment allowlist with SV junction k-mers when provided.
    // SV breakpoint sequences (e.g. deletion junctions, duplication junctions)
    // contain k-mers not present in the reference targets and would otherwise
    // be discarded by the allowlist filter.
    if let Some(ref junctions_path) = args.sv_junctions {
        let junctions = read_fasta(junctions_path)?;
        let junction_slices: Vec<&[u8]> =
            junctions.iter().map(|(_id, seq)| seq.as_slice()).collect();
        let junction_allowlist = build_allowlist(&junction_slices, k);
        let n_junction = junction_allowlist.len();
        allowlist.extend(junction_allowlist);
        eprintln!(
            "[index] sv_junctions: added {n_junction} junction k-mers ({} total)",
            allowlist.len()
        );
    }

    let n_target_kmers = allowlist.len() as u64;

    // ── 4. Build raw HashKmerIndex from molecules ─────────────────────────────
    let mut raw_index = HashKmerIndex::new();
    let reads = molecules_to_consensus_reads(&molecules);
    extract_all(&reads, k, &mut raw_index);

    // ── 5. Filter to allowlist entries and serialize ──────────────────────────
    let mut entries: Vec<KmerEntry> = raw_index
        .iter()
        .filter(|(&kmer, _)| allowlist.contains(&kmer))
        .map(|(&kmer, ev)| KmerEntry {
            kmer,
            n_molecules: ev.n_molecules,
            n_duplex: ev.n_duplex,
            n_simplex_fwd: ev.n_simplex_fwd,
            n_simplex_rev: ev.n_simplex_rev,
            min_base_error_prob: ev.min_base_error_prob,
            mean_base_error_prob: ev.mean_base_error_prob,
        })
        .collect();

    entries.sort_by_key(|e| e.kmer); // deterministic output order

    let n_kmers_observed = entries.len() as u64;
    let mean_molecule_depth = if n_kmers_observed > 0 {
        entries.iter().map(|e| e.n_molecules as f64).sum::<f64>() / n_kmers_observed as f64
    } else {
        0.0
    };

    write_bincode(&args.output, FileType::KmerIndex, &entries)?;

    // ── 6. Build and write QC JSON ────────────────────────────────────────────
    let qc = IndexQc {
        stage: "kmer_indexing".to_string(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        n_target_kmers,
        n_kmers_observed,
        mean_molecule_depth,
        passed: true,
    };

    let qc_path = args
        .output
        .parent()
        .unwrap_or(std::path::Path::new("."))
        .join("index_qc.json");
    write_qc(&qc_path, &qc)?;

    eprintln!(
        "[index] target_kmers={} observed={} mean_depth={:.2}",
        n_target_kmers, n_kmers_observed, mean_molecule_depth,
    );

    Ok(())
}

// ── Shared helpers ────────────────────────────────────────────────────────────

/// A parsed FASTA record as `(id, sequence)`.
pub type FastaRecord = (String, Vec<u8>);

/// Read a FASTA file and return `(id, sequence)` pairs.
///
/// Sequence bytes are uppercased on read.
pub fn read_fasta(path: &std::path::Path) -> Result<Vec<FastaRecord>, Box<dyn std::error::Error>> {
    let mut records = Vec::new();
    let mut reader = needletail::parse_fastx_file(path)?;
    while let Some(result) = reader.next() {
        let rec = result?;
        let id = String::from_utf8_lossy(rec.id())
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        let seq: Vec<u8> = rec.seq().iter().map(|&b| b.to_ascii_uppercase()).collect();
        records.push((id, seq));
    }
    Ok(records)
}

/// Convert assembled molecules to [`ConsensusReadInfo`] entries.
///
/// Prefers duplex consensus when available; falls back to forward or reverse
/// simplex consensus.
pub fn molecules_to_consensus_reads(molecules: &[Molecule]) -> Vec<ConsensusReadInfo> {
    let mut reads = Vec::new();
    for mol in molecules {
        if let Some(cr) = &mol.duplex_consensus {
            reads.push(ConsensusReadInfo {
                sequence: cr.sequence.clone(),
                per_base_error_prob: cr.per_base_error_prob.clone(),
                family_type: FamilyType::Duplex,
            });
        } else if let Some(cr) = &mol.consensus_fwd {
            reads.push(ConsensusReadInfo {
                sequence: cr.sequence.clone(),
                per_base_error_prob: cr.per_base_error_prob.clone(),
                family_type: FamilyType::from_family_size(cr.family_size),
            });
        } else if let Some(cr) = &mol.consensus_rev {
            reads.push(ConsensusReadInfo {
                sequence: cr.sequence.clone(),
                per_base_error_prob: cr.per_base_error_prob.clone(),
                family_type: FamilyType::from_family_size(cr.family_size),
            });
        }
    }
    reads
}

/// Reconstruct a [`HashKmerIndex`] from serialized [`KmerEntry`] records.
///
/// Used by downstream subcommands that read an index bincode file.
pub fn entries_to_hash_index(entries: &[KmerEntry]) -> HashKmerIndex {
    let mut index = HashKmerIndex::new();
    for e in entries {
        index.insert(
            e.kmer,
            MoleculeEvidence {
                n_molecules: e.n_molecules,
                n_duplex: e.n_duplex,
                n_simplex_fwd: e.n_simplex_fwd,
                n_simplex_rev: e.n_simplex_rev,
                min_base_error_prob: e.min_base_error_prob,
                mean_base_error_prob: e.mean_base_error_prob,
            },
        );
    }
    index
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as IoWrite;
    use std::path::PathBuf;

    use kam_assemble::assembler::{assemble_molecules, AssemblerConfig};
    use kam_assemble::parser::{parse_read_pair, ParseResult, ParserConfig};
    use kam_core::serialize::write_bincode;

    fn make_molecule() -> Molecule {
        // 5 UMI + 2 skip + 24 template = 31 bp
        let r1_seq = b"ACGTATGACGTACGTACGTACGTACGTACGT";
        let r2_seq = b"TGCATAGACGTACGTACGTACGTACGTACGT";
        let qual = vec![b'I'; r1_seq.len()];
        let config = ParserConfig::default();
        let pair =
            match parse_read_pair(r1_seq, &qual, r2_seq, &qual, &config).expect("parse error") {
                ParseResult::Ok(p) => p,
                ParseResult::Dropped { reason, detail } => {
                    panic!("unexpected drop: {reason:?}: {detail}");
                }
            };
        let (molecules, _) = assemble_molecules(vec![pair], &AssemblerConfig::default());
        molecules.into_iter().next().expect("at least one molecule")
    }

    fn write_fasta(path: &PathBuf, records: &[(&str, &str)]) {
        let mut f = std::fs::File::create(path).expect("create fasta");
        for (name, seq) in records {
            writeln!(f, ">{name}").unwrap();
            writeln!(f, "{seq}").unwrap();
        }
    }

    /// Integration test: molecules + targets FASTA → index.bin + index_qc.json.
    #[test]
    fn run_index_produces_output_files() {
        let dir = tempfile::tempdir().expect("tempdir");

        let mol = make_molecule();
        let molecules_path = dir.path().join("molecules.bin");
        write_bincode(&molecules_path, FileType::Molecules, &[mol]).expect("write molecules");

        let targets_path = dir.path().join("targets.fa");
        write_fasta(
            &targets_path,
            &[("target1", "ACGTACGTACGTACGTACGTACGTACGT")],
        );

        let output_path = dir.path().join("index.bin");

        let args = IndexArgs {
            input: molecules_path,
            targets: targets_path,
            output: output_path.clone(),
            kmer_size: 8,
            sv_junctions: None,
        };

        run_index(args).expect("run_index should succeed");

        assert!(output_path.exists(), "index.bin should exist");
        let qc_path = dir.path().join("index_qc.json");
        assert!(qc_path.exists(), "index_qc.json should exist");

        let qc_text = std::fs::read_to_string(&qc_path).unwrap();
        let _v: serde_json::Value =
            serde_json::from_str(&qc_text).expect("index_qc.json must be valid JSON");
    }

    /// Empty molecules file produces zero-entry index.
    #[test]
    fn run_index_empty_molecules_produces_empty_index() {
        let dir = tempfile::tempdir().expect("tempdir");

        let molecules_path = dir.path().join("molecules.bin");
        let empty: Vec<Molecule> = vec![];
        write_bincode(&molecules_path, FileType::Molecules, &empty).expect("write empty molecules");

        let targets_path = dir.path().join("targets.fa");
        write_fasta(&targets_path, &[("t1", "ACGTACGT")]);

        let output_path = dir.path().join("index.bin");

        let args = IndexArgs {
            input: molecules_path,
            targets: targets_path,
            output: output_path.clone(),
            kmer_size: 4,
            sv_junctions: None,
        };

        run_index(args).expect("run_index with empty molecules should succeed");
        assert!(output_path.exists());
    }
}
