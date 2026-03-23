//! Implementation of the `kam call` subcommand.
//!
//! Reads scored paths from a bincode file, calls variants against the reference
//! path per target, and writes output in the requested format(s).

use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;

use kam_call::caller::{call_variant, VariantFilter};
use kam_call::output::write_variants;
use kam_call::targeting::{
    apply_target_filter, apply_target_filter_with_tolerance, load_target_variants,
};
use kam_core::qc::{write_qc, CallQc};
use kam_core::serialize::read_bincode;

use crate::caller_config::caller_config_from_args;
use crate::cli::CallArgs;
use crate::commands::pathfind::ScoredPathRecord;
use crate::output::{format_extension, parse_output_formats};

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
    let caller_config = caller_config_from_args(&args);

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

                all_calls.push(call);
            }
        }
    }

    // ── 5b. Apply tumour-informed filter if --target-variants provided ─────────
    if let Some(ref vcf_path) = args.target_variants {
        let targets = load_target_variants(vcf_path)?;
        let tol = args.ti_position_tolerance as i64;
        if tol > 0 {
            apply_target_filter_with_tolerance(&mut all_calls, &targets, tol);
        } else {
            apply_target_filter(&mut all_calls, &targets);
        }
        eprintln!(
            "[call] tumour-informed filter applied: {} target variants loaded (position tolerance: {}bp)",
            targets.len(),
            tol,
        );
    }

    for call in &all_calls {
        if call.filter == VariantFilter::Pass {
            n_pass += 1;
        } else {
            n_filtered += 1;
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
        min_variant_specific_duplex: rec.min_variant_specific_duplex,
        min_simplex_fwd: rec.min_simplex_fwd,
        min_simplex_rev: rec.min_simplex_rev,
        mean_error_prob: rec.mean_error_prob,
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
            min_variant_specific_duplex: 0,
            min_simplex_fwd: n_molecules / 2,
            min_simplex_rev: n_molecules / 2,
            mean_error_prob: 0.001,
        }
    }

    fn default_call_args(
        paths_path: std::path::PathBuf,
        output_path: std::path::PathBuf,
    ) -> CallArgs {
        CallArgs {
            paths: paths_path,
            output: output_path,
            output_format: "tsv".to_string(),
            min_confidence: None,
            strand_bias_threshold: None,
            min_alt_molecules: None,
            min_alt_duplex: None,
            sv_min_confidence: None,
            sv_min_alt_molecules: None,
            max_vaf: None,
            target_variants: None,
            ti_position_tolerance: 0,
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
        let args = default_call_args(paths_path, output_path.clone());

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
        let args = default_call_args(paths_path, output_path.clone());

        run_call(args).expect("run_call with empty input should succeed");
        assert!(output_path.exists());
    }

    /// Multi-format output: requesting tsv,vcf produces both output files.
    #[test]
    fn run_call_multi_format_output() {
        let dir = tempfile::tempdir().expect("tempdir");

        let ref_rec = make_scored_path_record("TP53", b"ACGTACGT".to_vec(), true, 990);
        let alt_rec = make_scored_path_record("TP53", b"ACTTACGT".to_vec(), false, 10);

        let paths_path = dir.path().join("paths.bin");
        write_bincode(&paths_path, FileType::ScoredPaths, &[ref_rec, alt_rec])
            .expect("write paths");

        let output_path = dir.path().join("variants.tsv");

        let args = CallArgs {
            output_format: "tsv,vcf".to_string(),
            ..default_call_args(paths_path, output_path)
        };

        run_call(args).expect("run_call multi-format should succeed");

        assert!(
            dir.path().join("variants.tsv").exists(),
            "variants.tsv should exist"
        );
        assert!(
            dir.path().join("variants.vcf").exists(),
            "variants.vcf should exist"
        );
    }

    /// Tumour-informed filter: targeted calls pass, non-targeted calls are marked NotTargeted.
    ///
    /// Target IDs must use chrom:start-end format so the VCF matcher can extract
    /// coordinates. We use a generous position tolerance (1000 bp) to match without
    /// requiring exact sequence position alignment.
    #[test]
    fn run_call_target_filter_pass_and_not_targeted() {
        let dir = tempfile::tempdir().expect("tempdir");

        // Two targets: TP53 (will be in the tumour VCF) and BRCA1 (not in VCF).
        // Target IDs use chrom:start-end format required by parse_target_id.
        let ref_tp53 = make_scored_path_record("TP53:0-100", b"ACGTACGTACGT".to_vec(), true, 990);
        let alt_tp53 = make_scored_path_record("TP53:0-100", b"ACTTACGTACGT".to_vec(), false, 10);
        let ref_brca1 = make_scored_path_record("BRCA1:0-100", b"TTTTGGGGCCCC".to_vec(), true, 990);
        let alt_brca1 = make_scored_path_record("BRCA1:0-100", b"TTTTGGGGATCC".to_vec(), false, 10);

        let paths_path = dir.path().join("paths.bin");
        write_bincode(
            &paths_path,
            FileType::ScoredPaths,
            &[ref_tp53, alt_tp53, ref_brca1, alt_brca1],
        )
        .expect("write paths");

        // Write a minimal VCF with only the TP53 variant.
        let vcf_path = dir.path().join("targets.vcf");
        std::fs::write(
            &vcf_path,
            "##fileformat=VCFv4.2\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             TP53\t3\t.\tG\tT\t.\tPASS\t.\n",
        )
        .expect("write vcf");

        let output_path = dir.path().join("variants.vcf");

        let args = CallArgs {
            output_format: "vcf".to_string(),
            target_variants: Some(vcf_path),
            ti_position_tolerance: 1000,
            ..default_call_args(paths_path, output_path.clone())
        };

        run_call(args).expect("run_call with target filter should succeed");

        // Read the VCF output and verify filter assignments.
        let vcf_text = std::fs::read_to_string(&output_path).expect("read vcf");
        let data_lines: Vec<&str> = vcf_text.lines().filter(|l| !l.starts_with('#')).collect();

        // We must have two calls (one per target).
        assert_eq!(
            data_lines.len(),
            2,
            "expected 2 variant calls; output:\n{vcf_text}"
        );

        for line in &data_lines {
            let cols: Vec<&str> = line.split('\t').collect();
            let chrom = cols[0];
            let filter = cols[6];
            if chrom == "BRCA1" {
                assert_eq!(
                    filter, "NotTargeted",
                    "BRCA1 calls should be NotTargeted; output:\n{vcf_text}"
                );
            } else if chrom == "TP53" {
                assert_ne!(
                    filter, "NotTargeted",
                    "TP53 calls should not be NotTargeted; output:\n{vcf_text}"
                );
            }
        }
    }
}
