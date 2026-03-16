//! Implementation of the `kam assemble` subcommand.
//!
//! Reads paired FASTQ files, parses UMI/template structure, assembles
//! molecules, writes output in bincode format, and emits a QC JSON.

use kam_assemble::assembler::{assemble_molecules, AssemblerConfig};
use kam_assemble::consensus::ConsensusConfig;
use kam_assemble::io::read_fastq_pairs;
use kam_assemble::parser::ParserConfig;
use kam_core::qc::{write_qc, AssemblyQc};
use kam_core::serialize::{write_bincode, FileType};

use crate::cli::AssembleArgs;

/// Run the `assemble` subcommand end-to-end.
///
/// # Errors
///
/// Returns an error if any file I/O or serialization step fails.
pub fn run_assemble(args: AssembleArgs) -> Result<(), Box<dyn std::error::Error>> {
    // ── 1. Build parser config ────────────────────────────────────────────────
    let parser_config = ParserConfig {
        min_template_length: args.min_template_length.map(|v| v as usize),
        min_umi_quality: if args.min_umi_quality == 0 {
            None
        } else {
            Some(args.min_umi_quality)
        },
        ..ParserConfig::default()
    };

    // ── 2. Build assembler config ─────────────────────────────────────────────
    let assembler_config = AssemblerConfig {
        min_family_size: args.min_family_size as u8,
        consensus: ConsensusConfig::default(),
        ..AssemblerConfig::default()
    };

    // ── 3. Read FASTQ pairs ───────────────────────────────────────────────────
    let (read_pairs, parse_stats) =
        read_fastq_pairs(&args.r1, &args.r2, &parser_config)?;

    // ── 4. Assemble molecules ─────────────────────────────────────────────────
    let (molecules, mut assembly_stats) = assemble_molecules(read_pairs, &assembler_config);
    assembly_stats.parse_stats = parse_stats;

    // ── 5. Write molecules bincode ────────────────────────────────────────────
    write_bincode(&args.output, FileType::Molecules, &molecules)?;

    // ── 6. Build and write QC JSON ────────────────────────────────────────────
    let n_molecules = assembly_stats.n_molecules;
    let n_input = assembly_stats.parse_stats.n_processed;
    let n_dropped = n_input.saturating_sub(assembly_stats.parse_stats.n_passed);
    let duplex_fraction = if n_molecules > 0 {
        assembly_stats.n_duplex as f64 / n_molecules as f64
    } else {
        0.0
    };
    let mean_family_size = if n_molecules > 0 {
        (assembly_stats.parse_stats.n_passed as f64) / (n_molecules as f64)
    } else {
        0.0
    };

    let qc = AssemblyQc {
        stage: "molecule_assembly".to_string(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        n_input_read_pairs: n_input,
        n_molecules,
        n_duplex: assembly_stats.n_duplex,
        n_simplex_fwd: assembly_stats.n_simplex_fwd,
        n_simplex_rev: assembly_stats.n_simplex_rev,
        n_singletons: assembly_stats.n_singletons,
        duplex_fraction,
        mean_family_size,
        n_dropped_reads: n_dropped,
        passed: true,
    };

    // Derive QC path: same dir as output, named "assembly_qc.json".
    let qc_path = args
        .output
        .parent()
        .unwrap_or(std::path::Path::new("."))
        .join("assembly_qc.json");
    write_qc(&qc_path, &qc)?;

    // ── 7. Print summary to stderr ────────────────────────────────────────────
    eprintln!(
        "[assemble] input_pairs={} molecules={} duplex={} dropped={}",
        n_input, n_molecules, assembly_stats.n_duplex, n_dropped,
    );

    Ok(())
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as IoWrite;
    use std::path::PathBuf;

    /// Write a minimal FASTQ record set to `path`.
    fn write_fastq(path: &PathBuf, records: &[(&str, &str, &str)]) {
        let mut f = std::fs::File::create(path).expect("create fastq");
        for (name, seq, qual) in records {
            writeln!(f, "@{name}").unwrap();
            writeln!(f, "{seq}").unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{qual}").unwrap();
        }
    }

    /// Integration test: small FASTQ pair → bincode molecules + QC JSON.
    #[test]
    fn run_assemble_produces_output_files() {
        let dir = tempfile::tempdir().expect("tempdir");

        // Build read: 5 UMI + 2 skip + 20 template = 27 bp
        let r1_seq = "ACGTATGNNNNNNNNNNNNNNNNNNNN";
        let r2_seq = "TGCATAGNNNNNNNNNNNNNNNNNNNN";
        let qual = "I".repeat(r1_seq.len());

        let r1_path = dir.path().join("R1.fq");
        let r2_path = dir.path().join("R2.fq");
        write_fastq(&r1_path, &[("r1", r1_seq, &qual)]);
        write_fastq(&r2_path, &[("r1", r2_seq, &qual)]);

        let output_path = dir.path().join("molecules.bin");

        let args = AssembleArgs {
            r1: r1_path,
            r2: r2_path,
            output: output_path.clone(),
            chemistry: "twist-umi-duplex".to_string(),
            min_umi_quality: 0,
            min_family_size: 1,
            min_template_length: None,
            log_dir: None,
            log: vec![],
            threads: None,
        };

        run_assemble(args).expect("run_assemble should succeed");

        assert!(output_path.exists(), "molecules.bin should exist");
        let qc_path = dir.path().join("assembly_qc.json");
        assert!(qc_path.exists(), "assembly_qc.json should exist");

        // Verify QC JSON is valid JSON.
        let qc_text = std::fs::read_to_string(&qc_path).unwrap();
        let _v: serde_json::Value = serde_json::from_str(&qc_text)
            .expect("assembly_qc.json should be valid JSON");
    }

    /// Empty FASTQ input produces valid (empty) output.
    #[test]
    fn run_assemble_empty_input_produces_valid_output() {
        let dir = tempfile::tempdir().expect("tempdir");
        let r1_path = dir.path().join("R1.fq");
        let r2_path = dir.path().join("R2.fq");
        write_fastq(&r1_path, &[]);
        write_fastq(&r2_path, &[]);

        let output_path = dir.path().join("molecules.bin");

        let args = AssembleArgs {
            r1: r1_path,
            r2: r2_path,
            output: output_path.clone(),
            chemistry: "twist-umi-duplex".to_string(),
            min_umi_quality: 0,
            min_family_size: 1,
            min_template_length: None,
            log_dir: None,
            log: vec![],
            threads: None,
        };

        run_assemble(args).expect("run_assemble with empty input should succeed");
        assert!(output_path.exists());
    }
}
