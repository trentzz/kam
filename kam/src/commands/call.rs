//! Implementation of the `kam call` subcommand.
//!
//! Reads scored paths from a bincode file, calls variants against the reference
//! path per target, and writes output in the requested format(s).

use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;

use kam_call::caller::{call_variant, CallerConfig};
use kam_call::output::{write_variants, OutputFormat};
use kam_core::qc::{write_qc, CallQc};
use kam_core::serialize::read_bincode;

use crate::cli::CallArgs;
use crate::commands::pathfind::ScoredPathRecord;

/// Run the `call` subcommand end-to-end.
///
/// # Errors
///
/// Returns an error if file I/O or serialization fails.
pub fn run_call(args: CallArgs) -> Result<(), Box<dyn std::error::Error>> {
    // ── 1. Read scored paths ──────────────────────────────────────────────────
    let (_header, records): (_, Vec<ScoredPathRecord>) = read_bincode(&args.paths)?;

    // ── 2. Group paths by target ──────────────────────────────────────────────
    let mut by_target: HashMap<String, Vec<ScoredPathRecord>> = HashMap::new();
    for rec in records {
        by_target
            .entry(rec.target_id.clone())
            .or_default()
            .push(rec);
    }

    // ── 3. Build caller config from args ──────────────────────────────────────
    let caller_config = CallerConfig {
        min_confidence: args
            .min_confidence
            .unwrap_or(CallerConfig::default().min_confidence),
        strand_bias_threshold: args
            .strand_bias_threshold
            .unwrap_or(CallerConfig::default().strand_bias_threshold),
        ..CallerConfig::default()
    };

    // ── 4. Parse output formats ───────────────────────────────────────────────
    let formats = parse_output_formats(&args.output_format)?;

    // ── 5. Call variants per target ───────────────────────────────────────────
    let mut all_calls = Vec::new();
    let mut n_pass: u64 = 0;
    let mut n_filtered: u64 = 0;

    // Sort target IDs for deterministic output.
    let mut target_ids: Vec<String> = by_target.keys().cloned().collect();
    target_ids.sort();

    for target_id in &target_ids {
        let paths = &by_target[target_id];

        // Find the reference path (highest-evidence path matching the reference).
        let ref_path = paths.iter().find(|p| p.is_reference);

        if let Some(ref_rec) = ref_path {
            let ref_ev = record_to_path_evidence(ref_rec);

            // Call a variant for each non-reference path.
            for alt_rec in paths.iter().filter(|p| !p.is_reference) {
                let alt_ev = record_to_path_evidence(alt_rec);

                let call = call_variant(
                    target_id,
                    &ref_ev,
                    &alt_ev,
                    &ref_rec.sequence,
                    &alt_rec.sequence,
                    &caller_config,
                );

                use kam_call::caller::VariantFilter;
                if call.filter == VariantFilter::Pass {
                    n_pass += 1;
                } else {
                    n_filtered += 1;
                }

                all_calls.push(call);
            }
        }
    }

    let n_variants_called = all_calls.len() as u64;

    // ── 6. Write variant output file(s) ───────────────────────────────────────
    if formats.len() == 1 {
        let file = File::create(&args.output)?;
        let mut writer = BufWriter::new(file);
        write_variants(&all_calls, formats[0], &mut writer)?;
    } else {
        // Multiple formats: use the output path as base, add extension.
        let base = args.output.with_extension("");
        for &fmt in &formats {
            let ext = format_extension(fmt);
            let path = base.with_extension(ext);
            let file = File::create(&path)?;
            let mut writer = BufWriter::new(file);
            write_variants(&all_calls, fmt, &mut writer)?;
        }
    }

    // ── 7. Build and write QC JSON ────────────────────────────────────────────
    let qc = CallQc {
        stage: "variant_calling".to_string(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        n_variants_called,
        n_pass,
        n_filtered,
        passed: true,
    };

    let qc_path = args
        .output
        .parent()
        .unwrap_or(std::path::Path::new("."))
        .join("call_qc.json");
    write_qc(&qc_path, &qc)?;

    eprintln!(
        "[call] variants={} pass={} filtered={}",
        n_variants_called, n_pass, n_filtered,
    );

    Ok(())
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Convert a [`ScoredPathRecord`] to a [`PathEvidence`] for the caller.
fn record_to_path_evidence(rec: &ScoredPathRecord) -> kam_pathfind::score::PathEvidence {
    kam_pathfind::score::PathEvidence {
        min_molecules: rec.min_molecules,
        mean_molecules: rec.mean_molecules,
        min_duplex: rec.min_duplex,
        mean_duplex: rec.mean_duplex,
        total_simplex_fwd: rec.total_simplex_fwd,
        total_simplex_rev: rec.total_simplex_rev,
        mean_error_prob: rec.mean_error_prob,
    }
}

/// Parse a comma-separated list of format names into [`OutputFormat`] values.
fn parse_output_formats(s: &str) -> Result<Vec<OutputFormat>, Box<dyn std::error::Error>> {
    let mut formats = Vec::new();
    for token in s.split(',') {
        let fmt = match token.trim().to_ascii_lowercase().as_str() {
            "tsv" => OutputFormat::Tsv,
            "csv" => OutputFormat::Csv,
            "json" => OutputFormat::Json,
            "vcf" => OutputFormat::Vcf,
            other => {
                return Err(format!("unknown output format: '{other}'").into());
            }
        };
        formats.push(fmt);
    }
    if formats.is_empty() {
        formats.push(OutputFormat::Tsv);
    }
    Ok(formats)
}

/// Return the file extension for a given output format.
fn format_extension(fmt: OutputFormat) -> &'static str {
    match fmt {
        OutputFormat::Tsv => "tsv",
        OutputFormat::Csv => "csv",
        OutputFormat::Json => "json",
        OutputFormat::Vcf => "vcf",
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    use kam_core::serialize::{write_bincode, FileType};

    fn make_scored_path_record(
        target_id: &str,
        seq: Vec<u8>,
        is_reference: bool,
        n_molecules: u32,
    ) -> ScoredPathRecord {
        ScoredPathRecord {
            target_id: target_id.to_string(),
            sequence: seq,
            is_reference,
            min_molecules: n_molecules,
            mean_molecules: n_molecules as f32,
            min_duplex: 0,
            mean_duplex: 0.0,
            total_simplex_fwd: n_molecules / 2,
            total_simplex_rev: n_molecules / 2,
            mean_error_prob: 0.001,
        }
    }

    /// Integration test: scored paths bincode → variant TSV + call_qc.json.
    #[test]
    fn run_call_produces_output_files() {
        let dir = tempfile::tempdir().expect("tempdir");

        let ref_rec = make_scored_path_record("TP53", b"ACGTACGT".to_vec(), true, 990);
        let alt_rec = make_scored_path_record("TP53", b"ACTTACGT".to_vec(), false, 10);

        let paths_path = dir.path().join("paths.bin");
        write_bincode(&paths_path, FileType::ScoredPaths, &[ref_rec, alt_rec])
            .expect("write paths");

        let output_path = dir.path().join("variants.tsv");

        let args = CallArgs {
            paths: paths_path,
            output: output_path.clone(),
            output_format: "tsv".to_string(),
            min_confidence: None,
            strand_bias_threshold: None,
        };

        run_call(args).expect("run_call should succeed");

        assert!(output_path.exists(), "variants.tsv should exist");
        let qc_path = dir.path().join("call_qc.json");
        assert!(qc_path.exists(), "call_qc.json should exist");

        let qc_text = std::fs::read_to_string(&qc_path).unwrap();
        let _v: serde_json::Value =
            serde_json::from_str(&qc_text).expect("call_qc.json must be valid JSON");
    }

    /// Empty paths file produces no variants.
    #[test]
    fn run_call_empty_paths_produces_no_variants() {
        let dir = tempfile::tempdir().expect("tempdir");

        let paths_path = dir.path().join("paths.bin");
        let empty: Vec<ScoredPathRecord> = vec![];
        write_bincode(&paths_path, FileType::ScoredPaths, &empty).expect("write empty paths");

        let output_path = dir.path().join("variants.tsv");

        let args = CallArgs {
            paths: paths_path,
            output: output_path.clone(),
            output_format: "tsv".to_string(),
            min_confidence: None,
            strand_bias_threshold: None,
        };

        run_call(args).expect("run_call with empty input should succeed");
        assert!(output_path.exists());
    }
}
