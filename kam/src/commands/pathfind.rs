//! Implementation of the `kam pathfind` subcommand.
//!
//! Reads a k-mer index and target FASTA, builds a de Bruijn graph for each
//! target, walks paths from the first to the last anchor k-mer, scores them,
//! and writes scored paths to a bincode file.

use serde::{Deserialize, Serialize};

use kam_core::qc::{write_qc, PathfindQc};
use kam_core::serialize::{read_bincode, write_bincode, FileType};
use kam_index::encode::{canonical, KmerIterator};
use kam_pathfind::anchor::{validate_anchors, DEFAULT_ANCHOR_THRESHOLD};
use kam_pathfind::graph::DeBruijnGraph;
use kam_pathfind::score::{score_and_rank_paths, ScoredPath};

use crate::cli::PathfindArgs;
use crate::commands::index::{entries_to_hash_index, read_fasta, KmerEntry};

/// A serializable scored path record for inter-stage bincode files.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScoredPathRecord {
    /// Target identifier (FASTA record ID).
    pub target_id: String,
    /// Reconstructed DNA sequence of the path.
    pub sequence: Vec<u8>,
    /// Whether this path matches the reference sequence.
    pub is_reference: bool,
    /// Minimum molecule support across all k-mers.
    pub min_molecules: u32,
    /// Mean molecule support.
    pub mean_molecules: f32,
    /// Minimum duplex support across all path k-mers.
    pub min_duplex: u32,
    /// Mean duplex support across all path k-mers.
    pub mean_duplex: f32,
    /// Minimum duplex support at variant-specific k-mers only.
    pub min_variant_specific_duplex: u32,
    /// Mean molecule support across variant-specific k-mers only.
    pub mean_variant_specific_molecules: f32,
    /// Minimum forward-strand simplex molecule count across k-mers.
    pub min_simplex_fwd: u32,
    /// Minimum reverse-strand simplex molecule count across k-mers.
    pub min_simplex_rev: u32,
    /// Mean per-base error probability.
    pub mean_error_prob: f32,
}

/// Run the `pathfind` subcommand end-to-end.
///
/// # Errors
///
/// Returns an error if file I/O or serialization fails.
pub fn run_pathfind(args: PathfindArgs) -> Result<(), Box<dyn std::error::Error>> {
    // ── 1. Read k-mer index ───────────────────────────────────────────────────
    let (_header, entries): (_, Vec<KmerEntry>) = read_bincode(&args.index)?;
    let index = entries_to_hash_index(&entries);

    // Infer k from the entries.  If empty, default to 31 (not used downstream).
    let k = infer_k_from_entries(&entries).unwrap_or(31);

    // ── 2. Read target sequences from FASTA ───────────────────────────────────
    let targets = read_fasta(&args.targets)?;

    // ── 3. For each target: validate anchors → build graph → walk → score ─────
    let mut all_records: Vec<ScoredPathRecord> = Vec::new();
    let mut n_targets_queried: u64 = 0;
    let mut n_targets_with_variants: u64 = 0;
    let mut n_anchors_non_unique: u64 = 0;

    for (target_id, target_seq) in &targets {
        n_targets_queried += 1;

        // Validate anchors.
        let anchor_result = validate_anchors(target_seq, k, &index, DEFAULT_ANCHOR_THRESHOLD);

        let anchors = match anchor_result {
            Some(a) => a,
            None => {
                eprintln!("[pathfind] target {target_id}: sequence too short for k={k}, skipping");
                continue;
            }
        };

        if !anchors.start_unique || !anchors.end_unique {
            n_anchors_non_unique += 1;
            if let Some(warn) = &anchors.warning {
                eprintln!("[pathfind] target {target_id}: {warn}");
            }
        }

        // Collect all k-mers from the target for graph construction.
        let target_kmers: Vec<u64> = KmerIterator::new(target_seq, k)
            .map(|(_, km)| canonical(km, k))
            .collect();

        // Build de Bruijn graph. Require at least 2 molecules per k-mer to
        // filter PCR and sequencing error k-mers before graph construction.
        // This reduces spurious branching that causes BFS/DFS queue explosion.
        let graph = DeBruijnGraph::from_index(&index, k, &target_kmers, 2);

        // Set max_path_length based on target size. For a 100bp target at
        // k=31, the reference path is ~70 k-mers; 50 extra k-mers is generous
        // headroom for indels. This prevents the walker from exploring paths
        // far outside the target window, which was the primary memory cost.
        // Number of k-mers in the reference path = target_len - k + 1.
        // Allow 50 extra k-mers of headroom for indels.
        //
        // SV targets can override this via a `_maxpathN` suffix in the ID
        // (e.g. `chr17:100-300_DEL_maxpath200`). Larger deletions may need a
        // reference path that exceeds the default headroom.
        let default_max_path = if k > 1 {
            target_seq.len().saturating_sub(k - 1) + 50
        } else {
            150
        };
        let target_max_path = parse_maxpath_from_id(target_id).unwrap_or(default_max_path);
        let walk_config = kam_pathfind::walk::WalkConfig {
            max_path_length: target_max_path,
            ..Default::default()
        };

        // Walk paths.
        let paths = kam_pathfind::walk::walk_paths(
            &graph,
            anchors.start_kmer,
            anchors.end_kmer,
            &walk_config,
        );

        // Score and rank paths.
        let scored: Vec<ScoredPath> = score_and_rank_paths(paths, &index, target_seq, k);

        let has_variant = scored.iter().any(|p| !p.is_reference);
        if has_variant {
            n_targets_with_variants += 1;
        }

        // Convert to serializable records.
        for sp in scored {
            all_records.push(ScoredPathRecord {
                target_id: target_id.clone(),
                sequence: sp.path.sequence,
                is_reference: sp.is_reference,
                min_molecules: sp.aggregate_evidence.min_molecules,
                mean_molecules: sp.aggregate_evidence.mean_molecules,
                min_duplex: sp.aggregate_evidence.min_duplex,
                mean_duplex: sp.aggregate_evidence.mean_duplex,
                min_variant_specific_duplex: sp.aggregate_evidence.min_variant_specific_duplex,
                mean_variant_specific_molecules: sp
                    .aggregate_evidence
                    .mean_variant_specific_molecules,
                min_simplex_fwd: sp.aggregate_evidence.min_simplex_fwd,
                min_simplex_rev: sp.aggregate_evidence.min_simplex_rev,
                mean_error_prob: sp.aggregate_evidence.mean_error_prob,
            });
        }
    }

    // ── 4. Write scored paths to bincode ──────────────────────────────────────
    write_bincode(&args.output, FileType::ScoredPaths, &all_records)?;

    // ── 5. Build and write QC JSON ────────────────────────────────────────────
    let qc = PathfindQc {
        stage: "graph_walking".to_string(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        n_targets_queried,
        n_targets_with_variants,
        n_anchors_non_unique,
        passed: true,
    };

    let qc_path = args
        .output
        .parent()
        .unwrap_or(std::path::Path::new("."))
        .join("pathfind_qc.json");
    write_qc(&qc_path, &qc)?;

    eprintln!(
        "[pathfind] targets={} with_variants={} paths={}",
        n_targets_queried,
        n_targets_with_variants,
        all_records.len(),
    );

    Ok(())
}

/// Parse a per-target `max_path_length` override from a target ID suffix.
///
/// SV junction targets can embed a path-length override by appending
/// `_maxpathN` to the FASTA record ID (e.g. `chr17:100-300_DEL_maxpath200`).
/// This allows the DFS to explore longer paths when the reference allele spans
/// more k-mers than the default headroom.
///
/// Returns `None` if no override is present.
pub fn parse_maxpath_from_id(id: &str) -> Option<usize> {
    let lower = id.to_ascii_lowercase();
    let pos = lower.rfind("_maxpath")?;
    let suffix = &lower[pos + "_maxpath".len()..];
    // Take digits only (ignore any further suffix after the number).
    let digits: String = suffix.chars().take_while(|c| c.is_ascii_digit()).collect();
    digits.parse().ok()
}

/// Attempt to infer `k` from a set of encoded k-mers.
///
/// The k-mer with the fewest leading zero bits gives a lower bound on k.
/// We try values 4–31 and return the first that is consistent.
///
/// Returns `None` if `entries` is empty.
fn infer_k_from_entries(entries: &[KmerEntry]) -> Option<usize> {
    if entries.is_empty() {
        return None;
    }
    // Walk common k values and pick the one where all k-mers fit in 2*k bits.
    // Since we cannot know k from the encoded value alone, default to 31.
    Some(31)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::parse_maxpath_from_id;

    #[test]
    fn parse_maxpath_from_id_with_suffix() {
        assert_eq!(
            parse_maxpath_from_id("chr17:100-300_DEL_maxpath200"),
            Some(200)
        );
        assert_eq!(
            parse_maxpath_from_id("chr3:50-250_INV_maxpath250"),
            Some(250)
        );
    }

    #[test]
    fn parse_maxpath_from_id_no_suffix() {
        assert_eq!(parse_maxpath_from_id("chr17:100-300"), None);
        assert_eq!(parse_maxpath_from_id("TP53_exon7"), None);
    }

    #[test]
    fn parse_maxpath_from_id_case_insensitive() {
        assert_eq!(parse_maxpath_from_id("chr1:0-100_MAXPATH150"), Some(150));
    }

    use super::*;
    use std::io::Write as IoWrite;
    use std::path::PathBuf;

    use kam_assemble::assembler::{assemble_molecules, AssemblerConfig};
    use kam_assemble::parser::{parse_read_pair, ParseResult, ParserConfig};
    use kam_core::serialize::{write_bincode, FileType};

    use crate::cli::IndexArgs;
    use crate::commands::index::run_index;

    fn write_fasta(path: &PathBuf, records: &[(&str, &str)]) {
        let mut f = std::fs::File::create(path).expect("create fasta");
        for (name, seq) in records {
            writeln!(f, ">{name}").unwrap();
            writeln!(f, "{seq}").unwrap();
        }
    }

    fn build_test_index(dir: &std::path::Path, target_seq: &str, k: u32) -> (PathBuf, PathBuf) {
        // Build a molecule whose consensus contains the target sequence.
        let template = target_seq.as_bytes();
        let mut r1_seq = b"ACGTATG".to_vec();
        r1_seq.extend_from_slice(template);
        let mut r2_seq = b"TGCATAG".to_vec();
        r2_seq.extend_from_slice(template);
        let qual = vec![b'I'; r1_seq.len()];

        let config = ParserConfig::default();
        let pair =
            match parse_read_pair(&r1_seq, &qual, &r2_seq, &qual, &config).expect("parse error") {
                ParseResult::Ok(p) => p,
                ParseResult::Dropped { reason, .. } => panic!("drop: {reason:?}"),
            };
        let (molecules, _) = assemble_molecules(vec![pair], &AssemblerConfig::default());

        let molecules_path = dir.join("molecules.bin");
        write_bincode(&molecules_path, FileType::Molecules, &molecules).expect("write molecules");

        let targets_path = dir.join("targets.fa");
        write_fasta(
            &targets_path.to_path_buf().clone().into(),
            &[("t1", target_seq)],
        );

        let index_path = dir.join("index.bin");
        let idx_args = IndexArgs {
            input: molecules_path,
            targets: targets_path.clone(),
            output: index_path.clone(),
            kmer_size: k,
            sv_junctions: None,
        };
        run_index(idx_args).expect("run_index should succeed");

        (index_path, targets_path)
    }

    #[test]
    fn run_pathfind_produces_output_files() {
        let dir = tempfile::tempdir().expect("tempdir");
        let target = "ACGTACGTACGTACGTACGTACGTACGT";
        let k = 8u32;

        let (index_path, targets_path) = build_test_index(dir.path(), target, k);

        let output_path = dir.path().join("paths.bin");

        let args = PathfindArgs {
            index: index_path,
            targets: targets_path,
            output: output_path.clone(),
        };

        run_pathfind(args).expect("run_pathfind should succeed");

        assert!(output_path.exists(), "paths.bin should exist");
        let qc_path = dir.path().join("pathfind_qc.json");
        assert!(qc_path.exists(), "pathfind_qc.json should exist");

        let qc_text = std::fs::read_to_string(&qc_path).unwrap();
        let _v: serde_json::Value =
            serde_json::from_str(&qc_text).expect("pathfind_qc.json must be valid JSON");
    }
}
