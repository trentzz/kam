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
use kam_call::caller::{call_variant, VariantFilter};
use kam_call::output::write_variants;
use kam_call::targeting::{
    apply_target_filter, apply_target_filter_with_tolerance, load_target_variants,
};
use kam_core::kmer::KmerIndex;
use kam_core::qc::{write_qc, AssemblyQc, CallQc, IndexQc, PathfindQc};
use kam_index::allowlist::build_allowlist;
use kam_index::encode::{canonical, encode_kmer, reverse_complement, KmerIterator};
use kam_index::extract::ConsensusReadInfo;
use kam_index::HashKmerIndex;
use kam_pathfind::anchor::{validate_anchors, DEFAULT_ANCHOR_THRESHOLD};
use kam_pathfind::graph::DeBruijnGraph;
use kam_pathfind::score::{score_and_rank_paths, ScoredPath};
use kam_pathfind::walk::{find_alt_paths_from_reference, walk_paths_biased, GraphPath, WalkConfig};

use crate::caller_config::caller_config_from_args;
use crate::cli::RunArgs;
use crate::commands::index::{molecules_to_consensus_reads, read_fasta};
use crate::commands::pathfind::parse_maxpath_from_id;
use crate::output::{format_extension, parse_output_formats};

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
    let mut allowlist = build_allowlist(&target_slices, k);

    // Augment allowlist with SV junction k-mers when provided.
    // Keep a separate canonical k-mer set for use in the reference-guided alt
    // path search: only explore branches matching known junction k-mers.
    // The raw junction sequences are kept for synthetic alt path construction
    // when graph-based junction traversal fails (e.g. DUP at low VAF where
    // only J[0..2] are covered by reads).
    let mut junction_canonical_kmers: std::collections::HashSet<u64> =
        std::collections::HashSet::new();
    let mut sv_junctions_data: Option<Vec<(String, Vec<u8>)>> = None;
    if let Some(ref junctions_path) = args.sv_junctions {
        let junctions = read_fasta(junctions_path)?;
        let junction_slices: Vec<&[u8]> =
            junctions.iter().map(|(_id, seq)| seq.as_slice()).collect();
        let junction_allowlist = build_allowlist(&junction_slices, k);
        let n_junction = junction_allowlist.len();
        junction_canonical_kmers = junction_allowlist.clone();
        allowlist.extend(junction_allowlist);
        eprintln!(
            "[run/index] sv_junctions: added {n_junction} junction k-mers ({} total)",
            allowlist.len()
        );
        sv_junctions_data = Some(junctions);
    }

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
    let mut n_walk_no_paths: u64 = 0;
    let mut n_walk_ref_only: u64 = 0;
    let mut n_walk_alt_found: u64 = 0;
    let mut n_start_not_in_graph: u64 = 0;
    let mut n_end_not_in_graph: u64 = 0;
    let mut n_soft_anchor_recovered: u64 = 0;

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

        // Soft anchor search: if the exact start/end anchor k-mer is absent from
        // the graph, scan inward up to ANCHOR_WINDOW positions to find the
        // nearest in-graph k-mer.  Missing exact anchors (18% of targets at
        // 2M reads) arise when reads do not cover the first/last k=31 bases of
        // the target window due to short insert sizes or coverage imbalance.
        //
        // After walking with the inward anchors, each path sequence is padded
        // with the skipped prefix/suffix from the reference, so downstream
        // scoring and calling compare against the full target sequence unchanged.
        const ANCHOR_WINDOW: usize = 10;

        let (start_offset, start_raw) =
            match find_soft_anchor(target_seq, k, &graph, false, ANCHOR_WINDOW) {
                Some(v) => v,
                None => {
                    n_walk_no_paths += 1;
                    all_scored.push((target_id.clone(), vec![]));
                    continue;
                }
            };

        let (end_offset, end_raw) =
            match find_soft_anchor(target_seq, k, &graph, true, ANCHOR_WINDOW) {
                Some(v) => v,
                None => {
                    n_walk_no_paths += 1;
                    all_scored.push((target_id.clone(), vec![]));
                    continue;
                }
            };

        // Guard: start and end anchors must be distinct k-mers.  Identical
        // k-mers arise when the target is extremely short (target_len < k+1)
        // combined with aggressive soft anchoring from both ends.  Walk handles
        // start==end as a single-kmer path, but that is meaningless for variant
        // calling, so skip here.  Overlapping anchor WINDOWS are fine; what
        // matters is that the encoded k-mers differ.
        if start_raw == end_raw {
            n_walk_no_paths += 1;
            all_scored.push((target_id.clone(), vec![]));
            continue;
        }

        // Track exact anchor misses (offset > 0 means the exact anchor was absent).
        if start_offset > 0 {
            n_start_not_in_graph += 1;
            n_soft_anchor_recovered += 1;
        }
        if end_offset > 0 {
            n_end_not_in_graph += 1;
            n_soft_anchor_recovered += 1;
        }

        // Number of k-mers in the effective (inward-trimmed) reference path.
        // The initial DFS uses a small headroom (50) — it only needs to find the
        // reference path and nearby short-indel variants.
        // The reference-guided alt search uses a larger headroom when sv-junctions
        // are provided, to accommodate tandem duplications (which add up to ~200bp
        // to the alt path length).
        // SV targets can override the base limit via a `_maxpathN` suffix in the ID.
        let effective_len = target_seq.len() - start_offset - end_offset;
        let default_max_path = if k > 1 {
            effective_len.saturating_sub(k - 1) + 50
        } else {
            150
        };
        let target_max_path = parse_maxpath_from_id(target_id).unwrap_or(default_max_path);
        // Alt path search allows 300 extra k-mers when sv-junctions are active.
        let alt_max_path = if args.sv_junctions.is_some() {
            target_max_path + 250
        } else {
            target_max_path
        };
        let walk_config = WalkConfig {
            max_path_length: target_max_path,
            ..Default::default()
        };

        // Use evidence-biased DFS: sort successors by molecule count descending
        // so the high-evidence reference path is found before max_paths is
        // exhausted by low-evidence error-k-mer branches.
        let raw_paths = walk_paths_biased(&graph, start_raw, end_raw, &walk_config, |raw_km| {
            index
                .get(canonical(raw_km, k))
                .map(|e| e.n_molecules)
                .unwrap_or(0)
        });

        // Categorise walk outcome.
        if raw_paths.is_empty() {
            n_walk_no_paths += 1;
        }

        // Pad each path sequence with the reference prefix/suffix that was
        // excluded from the walk due to soft anchoring.  This allows
        // score_and_rank_paths to compare against the full target_seq.
        let prefix = &target_seq[..start_offset];
        let suffix = &target_seq[target_seq.len() - end_offset..];
        let paths: Vec<GraphPath> = if start_offset == 0 && end_offset == 0 {
            raw_paths
        } else {
            raw_paths
                .into_iter()
                .map(|mut p| {
                    let mut padded =
                        Vec::with_capacity(prefix.len() + p.sequence.len() + suffix.len());
                    padded.extend_from_slice(prefix);
                    padded.extend_from_slice(&p.sequence);
                    padded.extend_from_slice(suffix);
                    p.sequence = padded;
                    p
                })
                .collect()
        };

        // Reference-guided alt path search.
        //
        // Exhaustive DFS (even evidence-biased) is dominated by single-base
        // error bubbles: each sequencing error along the ~170-k-mer reference
        // path generates a complete error path, exhausting max_paths=100
        // before the DFS backtracks to an early SV branch point (e.g. at
        // position 20/170 for a deletion). This second pass walks the
        // reference path and explicitly searches for high-evidence non-
        // reference branches at each position.
        // Clone the reference path k-mers before moving `paths`.
        let ref_path_kmers: Option<Vec<u64>> = paths
            .iter()
            .find(|p| p.sequence == target_seq.as_slice())
            .map(|p| p.kmers.clone());
        let mut paths = paths;
        if let Some(ref_kmers) = ref_path_kmers {
            // Branch filter: when junction k-mers are available, only explore
            // branches that are known junction k-mers (prevents running DFS
            // from every error-bubble branch, which causes 5-minute runtimes).
            // Without junction k-mers, fall back to a generous evidence
            // threshold.
            let use_junction_filter = !junction_canonical_kmers.is_empty();

            let guided = find_alt_paths_from_reference(
                &graph,
                &ref_kmers,
                end_raw,
                alt_max_path,
                1,
                |raw_km| {
                    let can_km = canonical(raw_km, k);
                    if use_junction_filter {
                        // Require the junction k-mer to have at least 3
                        // supporting molecules to avoid launching sub-DFS from
                        // spurious accidental k-mer collisions between the
                        // junction allowlist and the target graph, which would
                        // cause an exhaustive search over the full graph.
                        junction_canonical_kmers.contains(&can_km)
                            && index
                                .get(can_km)
                                .map(|e| e.n_molecules >= 3)
                                .unwrap_or(false)
                    } else {
                        index
                            .get(can_km)
                            .map(|e| e.n_molecules >= 10)
                            .unwrap_or(false)
                    }
                },
                |raw_km| {
                    index
                        .get(canonical(raw_km, k))
                        .map(|e| e.n_molecules)
                        .unwrap_or(0)
                },
            );
            // Synthetic fallback for DUP junctions: when graph-based traversal
            // finds no alt path (junction chain J[3..29] absent from raw graph
            // at low VAF), reconstruct the alt sequence directly from the
            // target and junction sequences. This bypasses the de Bruijn graph
            // for the junction portion; the path is scored via the canonical index.
            // Add guided alt paths not already present (by sequence).
            let mut existing_seqs: std::collections::HashSet<Vec<u8>> =
                paths.iter().map(|p| p.sequence.clone()).collect();
            for alt_path in guided {
                if existing_seqs.insert(alt_path.sequence.clone()) {
                    paths.push(alt_path);
                }
            }

            // Always try synthetic DUP alt paths from the junction file.
            // The graph-based traversal cannot find DUP paths at low VAF because
            // the full junction k-mer chain is not in the raw graph; graph-derived
            // paths (SNVs, chimeric bridges) may already be in `paths` but are
            // unrelated to the DUP. The synthetic path is added unconditionally so
            // it is scored alongside any graph-derived paths.
            if let Some(ref jdata) = sv_junctions_data {
                for (_junc_id, junc_seq) in jdata {
                    if let Some(alt_path) = synthesize_dup_alt_path(target_seq, junc_seq, k) {
                        if existing_seqs.insert(alt_path.sequence.clone()) {
                            paths.push(alt_path);
                        }
                    }
                }
            }
        }

        // Score against the canonical evidence index; score_path canonicalizes
        // each raw path k-mer before the lookup.
        let scored = score_and_rank_paths(paths, &index, target_seq, k);
        let has_variant = scored.iter().any(|p| !p.is_reference);
        if has_variant {
            n_targets_with_variants += 1;
            n_walk_alt_found += 1;
        } else if !scored.is_empty() {
            n_walk_ref_only += 1;
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
        "[run/pathfind] targets={} with_variants={} no_paths={} ref_only={} alt_found={} start_missing={} end_missing={} soft_recovered={} time_ms={}",
        n_targets_queried,
        n_targets_with_variants,
        n_walk_no_paths,
        n_walk_ref_only,
        n_walk_alt_found,
        n_start_not_in_graph,
        n_end_not_in_graph,
        n_soft_anchor_recovered,
        t_pathfind.elapsed().as_millis()
    );

    // ── Stage 4: Call ─────────────────────────────────────────────────────────
    let t_call = std::time::Instant::now();
    let caller_config = caller_config_from_args(&args);

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
        let tol = args.ti_position_tolerance as i64;
        if tol > 0 {
            apply_target_filter_with_tolerance(&mut all_calls, &target_set, tol);
        } else {
            apply_target_filter(&mut all_calls, &target_set);
        }
        eprintln!(
            "[run/call] tumour-informed filter applied: {} target variants loaded (position tolerance: {}bp)",
            target_set.len(),
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

/// Find the nearest in-graph anchor k-mer within `window` positions of the
/// nominal anchor.
///
/// When `from_end` is `false`, scans start positions `0, 1, ..., window`
/// (inward from the start of the target).  When `from_end` is `true`, scans
/// end positions `target_len - k`, `target_len - k - 1`, ... (inward from the
/// end of the target).
///
/// Returns `Some((offset_from_nominal, kmer))` for the first in-graph k-mer
/// found, or `None` if no k-mer within the window is present in the graph.
fn find_soft_anchor(
    target_seq: &[u8],
    k: usize,
    graph: &DeBruijnGraph,
    from_end: bool,
    window: usize,
) -> Option<(usize, u64)> {
    let len = target_seq.len();
    for off in 0..=window {
        let pos = if from_end {
            len.checked_sub(k + off)?
        } else {
            off
        };
        if pos + k > len {
            break;
        }
        if let Some(km) = encode_kmer(&target_seq[pos..pos + k]) {
            if graph.contains(km) {
                return Some((off, km));
            }
        }
    }
    None
}

/// Construct a synthetic [`GraphPath`] for a tandem duplication alt allele.
///
/// At low VAF, reads barely cross the D1→D2 boundary, covering only the first
/// few junction k-mers (J[0..2]) but not J[3..29].  Graph-based traversal
/// therefore cannot find the alt path.  This function bypasses the graph by
/// reconstructing the expected alt sequence directly from the target and
/// junction sequences.
///
/// The DUP junction sequence is `last_half bp of D` + `first_half bp of D`,
/// where D = target[dup_start..dup_end] is the duplicated region.
/// - `junction[0..half]`  = target[dup_end-half..dup_end]  → locates dup_end
/// - `junction[half..2*half]` = target[dup_start..dup_start+half] → locates dup_start
///
/// In the VCF representation, REF = target[dup_start] (anchor base) and
/// ALT = target[dup_start] + D.  The copy is therefore inserted right after
/// the anchor base, not at the end of D.  Alt sequence:
///   target[0..dup_start+1] + D + target[dup_start+1..]
///
/// This matches the alt genome produced by a typical SV simulator: the inserted
/// copy immediately follows the anchor base, and the original reference continues
/// from target[dup_start+1] after the copy.
///
/// The returned path uses raw (forward) k-mer encodings and is scored by the
/// caller via the canonical index.
///
/// Returns `None` if the junction does not match this target (either half is
/// absent from the target sequence, or dup_start ≥ dup_end).
fn synthesize_dup_alt_path(target_seq: &[u8], junc_seq: &[u8], k: usize) -> Option<GraphPath> {
    if junc_seq.len() < 2 || !junc_seq.len().is_multiple_of(2) {
        return None;
    }
    let half = junc_seq.len() / 2;

    // left half = last `half` bp of D → dup_end = pos + half
    let dup_end = find_subseq(target_seq, &junc_seq[..half])? + half;
    // right half = first `half` bp of D → dup_start = pos (anchor base)
    let dup_start = find_subseq(target_seq, &junc_seq[half..])?;

    if dup_start >= dup_end || dup_end > target_seq.len() {
        return None;
    }

    // The VCF inserts the copy right after the anchor (dup_start).
    // Alt = target[0..dup_start+1] + D + target[dup_start+1..]
    let anchor_end = dup_start + 1;
    if anchor_end > target_seq.len() {
        return None;
    }
    let dup_len = dup_end - dup_start;
    let mut alt_seq = Vec::with_capacity(target_seq.len() + dup_len);
    alt_seq.extend_from_slice(&target_seq[..anchor_end]);
    alt_seq.extend_from_slice(&target_seq[dup_start..dup_end]);
    alt_seq.extend_from_slice(&target_seq[anchor_end..]);

    let kmers: Vec<u64> = KmerIterator::new(&alt_seq, k).map(|(_, km)| km).collect();
    let length = kmers.len();
    if length == 0 {
        return None;
    }

    Some(GraphPath {
        kmers,
        sequence: alt_seq,
        length,
    })
}

/// Find the first occurrence of `needle` in `haystack`.
///
/// Returns the start index, or `None` if absent.
fn find_subseq(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || needle.len() > haystack.len() {
        return None;
    }
    haystack.windows(needle.len()).position(|w| w == needle)
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
            min_alt_duplex: None,
            sv_min_confidence: None,
            sv_min_alt_molecules: None,
            sv_strand_bias_threshold: 1.0,
            max_vaf: None,
            sv_junctions: None,
            target_variants: None,
            ti_position_tolerance: 0,
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
            min_alt_duplex: None,
            sv_min_confidence: None,
            sv_min_alt_molecules: None,
            sv_strand_bias_threshold: 1.0,
            max_vaf: None,
            sv_junctions: None,
            target_variants: None,
            ti_position_tolerance: 0,
            output_format: "tsv".to_string(),
            qc_output: None,
            log_dir: None,
            log: vec![],
            threads: None,
        };

        run_pipeline(args).expect("run_pipeline with empty input should succeed");
        assert!(output_dir.join("assembly_qc.json").exists());
    }

    /// Multi-format output: requesting tsv,vcf produces both output files.
    #[test]
    fn run_pipeline_multi_format_output() {
        let dir = tempfile::tempdir().expect("tempdir");

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
            min_alt_duplex: None,
            sv_min_confidence: None,
            sv_min_alt_molecules: None,
            max_vaf: None,
            sv_junctions: None,
            target_variants: None,
            ti_position_tolerance: 0,
            sv_strand_bias_threshold: 1.0,
            output_format: "tsv,vcf".to_string(),
            qc_output: None,
            log_dir: None,
            log: vec![],
            threads: None,
        };

        run_pipeline(args).expect("run_pipeline multi-format should succeed");

        assert!(
            output_dir.join("variants.tsv").exists(),
            "variants.tsv should exist"
        );
        assert!(
            output_dir.join("variants.vcf").exists(),
            "variants.vcf should exist"
        );
    }
}
