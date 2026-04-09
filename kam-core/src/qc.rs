//! QC JSON output per pipeline stage.
//!
//! Each pipeline stage produces a structured QC struct that is serialized to
//! JSON for Nextflow validation between stages.
//!
//! ## Example
//!
//! ```no_run
//! use kam_core::qc::{AssemblyQc, write_qc};
//! use std::path::Path;
//!
//! let qc = AssemblyQc {
//!     stage: "molecule_assembly".to_string(),
//!     version: "0.1.0".to_string(),
//!     n_input_read_pairs: 1000,
//!     n_molecules: 500,
//!     n_duplex: 400,
//!     n_simplex_fwd: 50,
//!     n_simplex_rev: 30,
//!     n_singletons: 20,
//!     duplex_fraction: 0.8,
//!     mean_family_size: 2.0,
//!     n_dropped_reads: 10,
//!     passed: true,
//! };
//! write_qc(Path::new("/tmp/assembly_qc.json"), &qc).unwrap();
//! ```

use serde::Serialize;
use std::path::Path;

/// QC output from the assembly stage.
#[derive(Debug, Clone, Serialize)]
pub struct AssemblyQc {
    /// Stage identifier: `"molecule_assembly"`
    pub stage: String,
    /// Pipeline version string
    pub version: String,
    /// Total read pairs seen from the input FASTQ
    pub n_input_read_pairs: u64,
    /// Number of distinct molecules after UMI grouping
    pub n_molecules: u64,
    /// Molecules with duplex (both strands) support
    pub n_duplex: u64,
    /// Molecules with only forward-strand simplex support
    pub n_simplex_fwd: u64,
    /// Molecules with only reverse-strand simplex support
    pub n_simplex_rev: u64,
    /// Molecules that appeared as a single read pair (no family)
    pub n_singletons: u64,
    /// Fraction of molecules that are duplex (`n_duplex / n_molecules`)
    pub duplex_fraction: f64,
    /// Mean reads per molecule family
    pub mean_family_size: f64,
    /// Read pairs dropped before molecule construction (e.g. QC fail)
    pub n_dropped_reads: u64,
    /// Whether this stage passed its QC thresholds
    pub passed: bool,
}

/// QC output from the k-mer indexing stage.
#[derive(Debug, Clone, Serialize)]
pub struct IndexQc {
    /// Stage identifier: `"kmer_indexing"`
    pub stage: String,
    /// Pipeline version string
    pub version: String,
    /// Number of k-mers in the target allowlist
    pub n_target_kmers: u64,
    /// Number of distinct target k-mers actually observed in molecules
    pub n_kmers_observed: u64,
    /// Mean molecule depth across observed k-mers
    pub mean_molecule_depth: f64,
    /// Whether this stage passed its QC thresholds
    pub passed: bool,
}

/// QC output from the graph pathfinding stage.
#[derive(Debug, Clone, Serialize)]
pub struct PathfindQc {
    /// Stage identifier: `"graph_walking"`
    pub stage: String,
    /// Pipeline version string
    pub version: String,
    /// Number of target loci queried
    pub n_targets_queried: u64,
    /// Targets for which at least one variant path was found
    pub n_targets_with_variants: u64,
    /// Targets whose anchors were not unique in the graph
    pub n_anchors_non_unique: u64,
    /// Number of targets where the DFS walk hit the expansion budget and terminated early.
    pub n_targets_walk_budget_exceeded: u64,
    /// Whether this stage passed its QC thresholds
    pub passed: bool,
}

/// QC output from the variant calling stage.
#[derive(Debug, Clone, Serialize)]
pub struct CallQc {
    /// Stage identifier: `"variant_calling"`
    pub stage: String,
    /// Pipeline version string
    pub version: String,
    /// Total variants emitted before filtering
    pub n_variants_called: u64,
    /// Variants that passed all filters
    pub n_pass: u64,
    /// Variants that were filtered out
    pub n_filtered: u64,
    /// Whether this stage passed its QC thresholds
    pub passed: bool,
}

/// Write any QC struct to a JSON file at `path`.
///
/// The file is created (or truncated) and written atomically per OS
/// buffering rules. Pretty-printed JSON is used so files are human-readable.
///
/// # Errors
///
/// Returns an error if the file cannot be created or if serialization fails.
///
/// # Example
///
/// ```
/// use kam_core::qc::{CallQc, write_qc};
/// use std::path::Path;
///
/// let qc = CallQc {
///     stage: "variant_calling".to_string(),
///     version: "0.1.0".to_string(),
///     n_variants_called: 10,
///     n_pass: 8,
///     n_filtered: 2,
///     passed: true,
/// };
/// let dir = tempfile::tempdir().unwrap();
/// let path = dir.path().join("call_qc.json");
/// write_qc(&path, &qc).unwrap();
/// ```
pub fn write_qc<T: Serialize>(path: &Path, qc: &T) -> Result<(), Box<dyn std::error::Error>> {
    let file = std::fs::File::create(path)?;
    serde_json::to_writer_pretty(file, qc)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn assembly_qc_fixture() -> AssemblyQc {
        AssemblyQc {
            stage: "molecule_assembly".to_string(),
            version: "0.1.0".to_string(),
            n_input_read_pairs: 1000,
            n_molecules: 500,
            n_duplex: 400,
            n_simplex_fwd: 50,
            n_simplex_rev: 30,
            n_singletons: 20,
            duplex_fraction: 0.8,
            mean_family_size: 2.0,
            n_dropped_reads: 10,
            passed: true,
        }
    }

    /// write_qc produces a non-empty file containing valid JSON.
    #[test]
    fn test_write_qc_produces_valid_json() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("assembly_qc.json");
        let qc = assembly_qc_fixture();
        write_qc(&path, &qc).expect("write_qc should succeed");

        let contents = fs::read_to_string(&path).expect("file should be readable");
        assert!(!contents.is_empty(), "JSON output must not be empty");

        // Confirm it round-trips through serde_json
        let _value: serde_json::Value =
            serde_json::from_str(&contents).expect("output must be valid JSON");
    }

    /// AssemblyQc serializes with all expected fields present.
    #[test]
    fn test_assembly_qc_serializes_all_fields() {
        let qc = assembly_qc_fixture();
        let json = serde_json::to_string(&qc).expect("serialization must succeed");
        let value: serde_json::Value =
            serde_json::from_str(&json).expect("must deserialize to Value");

        assert!(value.get("stage").is_some(), "stage field must be present");
        assert!(value.get("version").is_some());
        assert!(value.get("n_input_read_pairs").is_some());
        assert!(value.get("n_molecules").is_some());
        assert!(value.get("n_duplex").is_some());
        assert!(value.get("n_simplex_fwd").is_some());
        assert!(value.get("n_simplex_rev").is_some());
        assert!(value.get("n_singletons").is_some());
        assert!(value.get("duplex_fraction").is_some());
        assert!(value.get("mean_family_size").is_some());
        assert!(value.get("n_dropped_reads").is_some());
        assert!(value.get("passed").is_some());
    }

    /// Deserialized JSON has the correct stage name.
    #[test]
    fn test_deserialized_stage_name() {
        let qc = assembly_qc_fixture();
        let json = serde_json::to_string(&qc).expect("serialization must succeed");
        let value: serde_json::Value =
            serde_json::from_str(&json).expect("must deserialize to Value");

        assert_eq!(
            value["stage"].as_str().unwrap(),
            "molecule_assembly",
            "stage must be 'molecule_assembly'"
        );
    }

    /// duplex_fraction is computed and serialized correctly.
    #[test]
    fn test_duplex_fraction_correct() {
        // 400 duplex out of 500 molecules = 0.8
        let n_molecules: u64 = 500;
        let n_duplex: u64 = 400;
        let duplex_fraction = n_duplex as f64 / n_molecules as f64;

        let qc = AssemblyQc {
            stage: "molecule_assembly".to_string(),
            version: "0.1.0".to_string(),
            n_input_read_pairs: 1000,
            n_molecules,
            n_duplex,
            n_simplex_fwd: 50,
            n_simplex_rev: 30,
            n_singletons: 20,
            duplex_fraction,
            mean_family_size: 2.0,
            n_dropped_reads: 10,
            passed: true,
        };

        let json = serde_json::to_string(&qc).expect("serialization must succeed");
        let value: serde_json::Value =
            serde_json::from_str(&json).expect("must deserialize to Value");

        let serialized_fraction = value["duplex_fraction"]
            .as_f64()
            .expect("duplex_fraction must be a float");
        assert!(
            (serialized_fraction - 0.8).abs() < 1e-9,
            "duplex_fraction should be 0.8, got {serialized_fraction}"
        );
    }
}
