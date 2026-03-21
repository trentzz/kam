//! Implementation of the `kam run` full-pipeline subcommand.
//!
//! Runs all stages — assemble → index → pathfind → call — in a single process
//! with zero-copy data passing between stages (no bincode serialization
//! between stages).

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::BufWriter;

use kam_assemble::assembler::{assemble_molecules, AssemblerConfig};
use kam_assemble::consensus::ConsensusConfig;
use kam_assemble::io::read_fastq_pairs;
use kam_assemble::parser::ParserConfig;
use kam_call::caller::{call_variant, CallerConfig, VariantFilter};
use kam_call::output::{write_variants, OutputFormat};
use kam_call::targeting::{apply_target_filter, load_target_variants};
use kam_core::kmer::KmerIndex;
use kam_core::qc::{write_qc, AssemblyQc, CallQc, IndexQc, PathfindQc};
use kam_index::allowlist::build_allowlist;
use kam_index::encode::{canonical, encode_kmer, reverse_complement, KmerIterator};
use kam_index::extract::ConsensusReadInfo;
use kam_index::HashKmerIndex;
use kam_pathfind::anchor::{validate_anchors, DEFAULT_ANCHOR_THRESHOLD};
use kam_pathfind::graph::DeBruijnGraph;
use kam_pathfind::score::{score_and_rank_paths, ScoredPath};
use kam_pathfind::walk::WalkConfig;

use crate::cli::RunArgs;
use crate::commands::index::{molecules_to_consensus_reads, read_fasta};

/// Run the full pipeline end-to-end in memory (zero-copy hot path).
///
/// # Errors
///
/// Returns an error if any file I/O step fails.
pub fn run_pipeline(args: RunArgs) -> Result<(), Box<dyn std::error::Error>> {
    fs::create_dir_all(&args.output_dir)?;

    let k = args.kmer_size as usize;
    let t_total = std::time::Instant::now();

    // ── Stage 1: Assemble ─────────────────────────────────────────────────────
    let t_assemble = std::time::Instant::now();
    let parser_config = ParserConfig {
        min_template_length: args.min_template_length.map(|v| v as usize),
        min_umi_quality: if args.min_umi_quality == 0 {
            None
        } else {
            Some(args.min_umi_quality)
        },
        ..ParserConfig::default()
    };

    let assembler_config = AssemblerConfig {
        min_family_size: args.min_family_size as u8,
        consensus: ConsensusConfig::default(),
        ..AssemblerConfig::default()
    };

    let (read_pairs, parse_stats) = read_fastq_pairs(&args.r1, &args.r2, &parser_config)?;
    let (molecules, mut assembly_stats) = assemble_molecules(read_pairs, &assembler_config);
    assembly_stats.parse_stats = parse_stats;

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

    let assembly_qc = AssemblyQc {
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
    write_qc(&args.output_dir.join("assembly_qc.json"), &assembly_qc)?;
    eprintln!(
        "[run/assemble] molecules={} duplex={} time_ms={}",
        n_molecules,
        assembly_stats.n_duplex,
        t_assemble.elapsed().as_millis()
    );

    // ── Stage 2: Index ────────────────────────────────────────────────────────
    let t_index = std::time::Instant::now();
    let targets = read_fasta(&args.targets)?;
    let target_slices: Vec<&[u8]> = targets.iter().map(|(_id, seq)| seq.as_slice()).collect();
    let allowlist = build_allowlist(&target_slices, k);
    let n_target_kmers = allowlist.len() as u64;

    // Two-pass indexing: include ALL k-mers from molecules that overlap a
    // target region, not just exact target k-mers.
    //
    // Pass 1 (allowlist check): a molecule "overlaps a target" if ANY of its
    //   k-mers is in the allowlist.  This identifies on-target molecules.
    //
    // Pass 2 (full extraction): index every k-mer from on-target molecules,
    //   including variant k-mers that differ from the reference.  These
    //   non-allowlist k-mers are essential for de Bruijn graph variant paths.
    //
    // Off-target molecules (no k-mer in allowlist) are skipped for efficiency.
    let reads: Vec<ConsensusReadInfo> = molecules_to_consensus_reads(&molecules);
    let mut index = HashKmerIndex::new();
    for read in &reads {
        let overlaps_target = KmerIterator::new(&read.sequence, k)
            .any(|(_, km)| allowlist.contains(&canonical(km, k)));
        if overlaps_target {
            kam_index::extract::extract_and_index(read, k, &mut index);
        }
    }

    let n_kmers_observed = index.len() as u64;
    let mean_molecule_depth = if n_kmers_observed > 0 {
        index
            .iter()
            .map(|(_, ev)| ev.n_molecules as f64)
            .sum::<f64>()
            / n_kmers_observed as f64
    } else {
        0.0
    };

    // Build a raw (non-canonical) k-mer index for de Bruijn graph construction.
    //
    // The canonical index stores min(kmer, rev_comp(kmer)) for evidence
    // aggregation.  This is correct for counting but breaks the de Bruijn
    // graph: two consecutive k-mers from a sequence share a (k-1)-base overlap
    // only when both are in their original orientation.  After canonicalization
    // one or both can be flipped to their reverse complement, destroying the
    // suffix/prefix relationship and preventing edge formation.
    //
    // The raw index stores k-mers in the original orientation AND their reverse
    // complements so that every read direction is represented.  The graph is
    // then built from this raw index; scoring still uses the canonical index.
    // raw_graph_index is a presence map: stores raw (non-canonical) k-mers with
    // n_molecules=1 as a presence indicator.  HashKmerIndex::insert overwrites on
    // duplicate, so repeated insertions leave n_molecules=1.  This index is used
    // only for de Bruijn graph topology — evidence comes from the canonical index.
    let mut raw_graph_index = HashKmerIndex::new();
    let present = kam_core::kmer::MoleculeEvidence {
        n_molecules: 1,
        ..Default::default()
    };
    for read in &reads {
        let overlaps_target = KmerIterator::new(&read.sequence, k)
            .any(|(_, km)| allowlist.contains(&canonical(km, k)));
        if overlaps_target {
            for (_, raw_km) in KmerIterator::new(&read.sequence, k) {
                raw_graph_index.insert(raw_km, present.clone());
                raw_graph_index.insert(reverse_complement(raw_km, k), present.clone());
            }
        }
    }

    let index_qc = IndexQc {
        stage: "kmer_indexing".to_string(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        n_target_kmers,
        n_kmers_observed,
        mean_molecule_depth,
        passed: true,
    };
    write_qc(&args.output_dir.join("index_qc.json"), &index_qc)?;
    eprintln!(
        "[run/index] target_kmers={} observed={} time_ms={}",
        n_target_kmers,
        n_kmers_observed,
        t_index.elapsed().as_millis()
    );

    // ── Stage 3: Pathfind ─────────────────────────────────────────────────────
    let t_pathfind = std::time::Instant::now();

    // Build a single de Bruijn graph from all on-target raw k-mers.
    // The graph must include variant k-mers (which differ from the reference),
    // so it is built from ALL k-mers in raw_graph_index, not just reference k-mers.
    // min_molecules=1 accepts every k-mer observed in at least one molecule.
    // The walk is bounded per-target by max_path_length to prevent cross-target paths.
    let all_raw_kmers: Vec<u64> = raw_graph_index.iter().map(|(&km, _)| km).collect();
    let graph = DeBruijnGraph::from_index(&raw_graph_index, k, &all_raw_kmers, 1);

    let mut all_scored: Vec<(String, Vec<ScoredPath>)> = Vec::new();
    let mut n_targets_queried: u64 = 0;
    let mut n_targets_with_variants: u64 = 0;
    let mut n_anchors_non_unique: u64 = 0;

    for (target_id, target_seq) in &targets {
        n_targets_queried += 1;

        // validate_anchors checks canonical anchor uniqueness in the canonical
        // evidence index — used to warn about repeat regions.
        let anchor_result = validate_anchors(target_seq, k, &index, DEFAULT_ANCHOR_THRESHOLD);

        let anchors = match anchor_result {
            Some(a) => a,
            None => {
                eprintln!("[run/pathfind] target {target_id}: too short for k={k}, skipping");
                continue;
            }
        };

        if !anchors.start_unique || !anchors.end_unique {
            n_anchors_non_unique += 1;
        }

        // Use raw (non-canonical) anchors from the forward strand of the target.
        // The graph is built from raw k-mers so the walk starts and ends on
        // raw k-mer nodes.
        let start_raw = match encode_kmer(&target_seq[..k]) {
            Some(km) => km,
            None => {
                eprintln!("[run/pathfind] target {target_id}: start anchor contains N, skipping");
                continue;
            }
        };
        let end_raw = match encode_kmer(&target_seq[target_seq.len() - k..]) {
            Some(km) => km,
            None => {
                eprintln!("[run/pathfind] target {target_id}: end anchor contains N, skipping");
                continue;
            }
        };

        // Number of k-mers in the reference path = target_len - k + 1.
        // Allow 50 extra k-mers of headroom for indels.
        let target_max_path = if k > 1 {
            target_seq.len().saturating_sub(k - 1) + 50
        } else {
            150
        };
        let walk_config = WalkConfig {
            max_path_length: target_max_path,
            ..Default::default()
        };

        let paths = kam_pathfind::walk::walk_paths(&graph, start_raw, end_raw, &walk_config);

        // Score against the canonical evidence index; score_path canonicalizes
        // each raw path k-mer before the lookup.
        let scored = score_and_rank_paths(paths, &index, target_seq, k);
        let has_variant = scored.iter().any(|p| !p.is_reference);
        if has_variant {
            n_targets_with_variants += 1;
        }

        all_scored.push((target_id.clone(), scored));
    }

    let pathfind_qc = PathfindQc {
        stage: "graph_walking".to_string(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        n_targets_queried,
        n_targets_with_variants,
        n_anchors_non_unique,
        passed: true,
    };
    write_qc(&args.output_dir.join("pathfind_qc.json"), &pathfind_qc)?;
    eprintln!(
        "[run/pathfind] targets={} with_variants={} time_ms={}",
        n_targets_queried,
        n_targets_with_variants,
        t_pathfind.elapsed().as_millis()
    );

    // ── Stage 4: Call ─────────────────────────────────────────────────────────
    let t_call = std::time::Instant::now();
    let caller_config = CallerConfig {
        min_confidence: args
            .min_confidence
            .unwrap_or(CallerConfig::default().min_confidence),
        strand_bias_threshold: args
            .strand_bias_threshold
            .unwrap_or(CallerConfig::default().strand_bias_threshold),
        min_alt_molecules: args
            .min_alt_molecules
            .unwrap_or(CallerConfig::default().min_alt_molecules),
        max_vaf: args.max_vaf.or(CallerConfig::default().max_vaf),
        ..CallerConfig::default()
    };

    let mut all_calls = Vec::new();
    let mut n_pass: u64 = 0;
    let mut n_filtered: u64 = 0;

    // Build a map from target_id → target_seq for reference lookup.
    let target_map: HashMap<&str, &[u8]> = targets
        .iter()
        .map(|(id, seq)| (id.as_str(), seq.as_slice()))
        .collect();

    for (target_id, scored_paths) in &all_scored {
        let ref_path = scored_paths.iter().find(|p| p.is_reference);
        let ref_seq = target_map.get(target_id.as_str()).copied().unwrap_or(&[]);

        if let Some(ref_sp) = ref_path {
            for alt_sp in scored_paths.iter().filter(|p| !p.is_reference) {
                let call = call_variant(
                    target_id,
                    &ref_sp.aggregate_evidence,
                    &alt_sp.aggregate_evidence,
                    ref_seq,
                    &alt_sp.path.sequence,
                    &caller_config,
                );
                all_calls.push(call);
            }
        }
    }

    // Apply tumour-informed filter if --target-variants provided.
    if let Some(ref vcf_path) = args.target_variants {
        let target_set = load_target_variants(vcf_path)?;
        apply_target_filter(&mut all_calls, &target_set);
        eprintln!(
            "[run/call] tumour-informed filter applied: {} target variants loaded",
            target_set.len()
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

    let call_qc = CallQc {
        stage: "variant_calling".to_string(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        n_variants_called,
        n_pass,
        n_filtered,
        passed: true,
    };
    write_qc(&args.output_dir.join("call_qc.json"), &call_qc)?;
    eprintln!(
        "[run/call] variants={} pass={} filtered={} time_ms={}",
        n_variants_called,
        n_pass,
        n_filtered,
        t_call.elapsed().as_millis()
    );

    // ── Write final variant output ────────────────────────────────────────────
    let t_output = std::time::Instant::now();
    let formats = parse_output_formats(&args.output_format)?;
    let base_path = args.output_dir.join("variants");

    if formats.len() == 1 {
        let ext = format_extension(formats[0]);
        let path = base_path.with_extension(ext);
        let file = File::create(&path)?;
        let mut writer = BufWriter::new(file);
        write_variants(&all_calls, formats[0], &mut writer)?;
    } else {
        for &fmt in &formats {
            let ext = format_extension(fmt);
            let path = base_path.with_extension(ext);
            let file = File::create(&path)?;
            let mut writer = BufWriter::new(file);
            write_variants(&all_calls, fmt, &mut writer)?;
        }
    }

    eprintln!(
        "[run] output time_ms={} total_ms={}",
        t_output.elapsed().as_millis(),
        t_total.elapsed().as_millis()
    );
    eprintln!(
        "[run] pipeline complete — output in {}",
        args.output_dir.display()
    );

    Ok(())
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Parse a comma-separated format string into [`OutputFormat`] values.
fn parse_output_formats(s: &str) -> Result<Vec<OutputFormat>, Box<dyn std::error::Error>> {
    let mut formats = Vec::new();
    for token in s.split(',') {
        let fmt = match token.trim().to_ascii_lowercase().as_str() {
            "tsv" => OutputFormat::Tsv,
            "csv" => OutputFormat::Csv,
            "json" => OutputFormat::Json,
            "vcf" => OutputFormat::Vcf,
            other => return Err(format!("unknown output format: '{other}'").into()),
        };
        formats.push(fmt);
    }
    if formats.is_empty() {
        formats.push(OutputFormat::Tsv);
    }
    Ok(formats)
}

/// Return the file extension string for an output format.
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
    use std::io::Write as IoWrite;
    use std::path::PathBuf;

    fn write_fastq(path: &PathBuf, records: &[(&str, &str, &str)]) {
        let mut f = std::fs::File::create(path).expect("create fastq");
        for (name, seq, qual) in records {
            writeln!(f, "@{name}").unwrap();
            writeln!(f, "{seq}").unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{qual}").unwrap();
        }
    }

    fn write_fasta(path: &PathBuf, records: &[(&str, &str)]) {
        let mut f = std::fs::File::create(path).expect("create fasta");
        for (name, seq) in records {
            writeln!(f, ">{name}").unwrap();
            writeln!(f, "{seq}").unwrap();
        }
    }

    /// Integration test: full pipeline with small synthetic data.
    #[test]
    fn run_pipeline_produces_output_files() {
        let dir = tempfile::tempdir().expect("tempdir");

        // Build read: 5 UMI + 2 skip + template = 31 bp
        let template = "ACGTACGTACGTACGTACGTACGT";
        let r1_seq = format!("ACGTATG{template}");
        let r2_seq = format!("TGCATAG{template}");
        let qual = "I".repeat(r1_seq.len());

        let r1_path = dir.path().join("R1.fq");
        let r2_path = dir.path().join("R2.fq");
        write_fastq(&r1_path, &[("read1", &r1_seq, &qual)]);
        write_fastq(&r2_path, &[("read1", &r2_seq, &qual)]);

        let targets_path = dir.path().join("targets.fa");
        write_fasta(&targets_path, &[("target1", template)]);

        let output_dir = dir.path().join("results");

        let args = RunArgs {
            r1: r1_path,
            r2: r2_path,
            targets: targets_path,
            output_dir: output_dir.clone(),
            chemistry: "twist-umi-duplex".to_string(),
            min_umi_quality: 0,
            min_family_size: 1,
            min_template_length: None,
            kmer_size: 8,
            min_confidence: None,
            strand_bias_threshold: None,
            min_alt_molecules: None,
            max_vaf: None,
            target_variants: None,
            output_format: "tsv".to_string(),
            qc_output: None,
            log_dir: None,
            log: vec![],
            threads: None,
        };

        run_pipeline(args).expect("run_pipeline should succeed");

        // Check all QC files exist.
        assert!(
            output_dir.join("assembly_qc.json").exists(),
            "assembly_qc.json"
        );
        assert!(output_dir.join("index_qc.json").exists(), "index_qc.json");
        assert!(
            output_dir.join("pathfind_qc.json").exists(),
            "pathfind_qc.json"
        );
        assert!(output_dir.join("call_qc.json").exists(), "call_qc.json");
        assert!(output_dir.join("variants.tsv").exists(), "variants.tsv");

        // Verify QC JSONs are valid.
        for name in &[
            "assembly_qc.json",
            "index_qc.json",
            "pathfind_qc.json",
            "call_qc.json",
        ] {
            let text = std::fs::read_to_string(output_dir.join(name)).unwrap();
            let _v: serde_json::Value = serde_json::from_str(&text)
                .unwrap_or_else(|e| panic!("{name} is not valid JSON: {e}"));
        }
    }

    /// Empty FASTQ input completes without error.
    #[test]
    fn run_pipeline_empty_input_completes() {
        let dir = tempfile::tempdir().expect("tempdir");

        let r1_path = dir.path().join("R1.fq");
        let r2_path = dir.path().join("R2.fq");
        write_fastq(&r1_path, &[]);
        write_fastq(&r2_path, &[]);

        let targets_path = dir.path().join("targets.fa");
        write_fasta(&targets_path, &[("t1", "ACGTACGT")]);

        let output_dir = dir.path().join("out");

        let args = RunArgs {
            r1: r1_path,
            r2: r2_path,
            targets: targets_path,
            output_dir: output_dir.clone(),
            chemistry: "twist-umi-duplex".to_string(),
            min_umi_quality: 0,
            min_family_size: 1,
            min_template_length: None,
            kmer_size: 4,
            min_confidence: None,
            strand_bias_threshold: None,
            min_alt_molecules: None,
            max_vaf: None,
            target_variants: None,
            output_format: "tsv".to_string(),
            qc_output: None,
            log_dir: None,
            log: vec![],
            threads: None,
        };

        run_pipeline(args).expect("run_pipeline with empty input should succeed");
        assert!(output_dir.join("assembly_qc.json").exists());
    }
}
