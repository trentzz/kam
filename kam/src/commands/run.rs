//! Implementation of the `kam run` full-pipeline subcommand.
//!
//! Runs all stages — assemble → index → pathfind → call — in a single process
//! with zero-copy data passing between stages (no bincode serialisation
//! between stages).
//!
//! Configuration priority: CLI flags > config file (`--config`) > built-in defaults.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::BufWriter;

use kam_assemble::assembler::assemble_molecules;
use kam_assemble::io::read_fastq_pairs;
use kam_call::caller::{call_variant, VariantFilter};
use kam_call::fusion::{call_fusion, parse_fusion_targets, FusionContext};
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

use crate::cli::RunArgs;
use crate::commands::index::{molecules_to_consensus_reads, read_fasta};
use crate::commands::pathfind::parse_maxpath_from_id;
use crate::config::KamConfig;
use crate::output::{format_extension, parse_output_formats};

/// Run the full pipeline end-to-end in memory (zero-copy hot path).
///
/// When `args.config` is provided, pipeline parameters are loaded from the
/// TOML file first, then any CLI flags override individual values. When
/// `args.config` is absent, all parameters come directly from the CLI.
///
/// # Errors
///
/// Returns an error if any file I/O step fails or required parameters are absent.
pub fn run_pipeline(args: RunArgs) -> Result<(), Box<dyn std::error::Error>> {
    // Build the unified config: load from file if --config given, then merge CLI.
    let cfg = build_config(&args)?;

    // Unwrap required paths — validate() already checked them.
    let r1 = cfg.input.r1.as_ref().expect("r1 validated");
    let r2 = cfg.input.r2.as_ref().expect("r2 validated");
    let targets_path = cfg.input.targets.as_ref().expect("targets validated");
    let output_dir = cfg
        .output
        .output_dir
        .as_ref()
        .expect("output_dir validated");

    fs::create_dir_all(output_dir)?;

    let k = cfg.kmer_size();
    let t_total = std::time::Instant::now();

    // ── Stage 1: Assemble ─────────────────────────────────────────────────────
    let t_assemble = std::time::Instant::now();
    let parser_config = cfg.to_parser_config();
    let assembler_config = cfg.to_assembler_config();

    let (read_pairs, parse_stats) = read_fastq_pairs(r1, r2, &parser_config)?;
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
    write_qc(&output_dir.join("assembly_qc.json"), &assembly_qc)?;
    eprintln!(
        "[run/assemble] molecules={} duplex={} time_ms={}",
        n_molecules,
        assembly_stats.n_duplex,
        t_assemble.elapsed().as_millis()
    );

    // ── Stage 2: Index ────────────────────────────────────────────────────────
    let t_index = std::time::Instant::now();
    let targets = read_fasta(targets_path)?;
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
    if let Some(ref junctions_path) = cfg.input.sv_junctions {
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

    // Augment allowlist with fusion target k-mers when provided.
    // Fusion targets are synthetic sequences; their k-mers must be in the
    // allowlist so that fusion-spanning reads are captured in the two-pass index.
    // Normal targets and fusion targets are kept separate: fusion targets are
    // processed after normal target calling so that partner depths are available.
    let fusion_targets_data = if let Some(ref fusion_path) = cfg.input.fusion_targets {
        let fts = parse_fusion_targets(fusion_path)?;
        let fusion_slices: Vec<&[u8]> = fts.iter().map(|t| t.sequence.as_slice()).collect();
        let fusion_allowlist = build_allowlist(&fusion_slices, k);
        let n_fusion = fusion_allowlist.len();
        allowlist.extend(fusion_allowlist);
        eprintln!(
            "[run/index] fusion_targets: loaded {} targets, added {n_fusion} k-mers ({} total)",
            fts.len(),
            allowlist.len()
        );
        fts
    } else {
        Vec::new()
    };

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
    //
    // Pass 1: count how many on-target molecules each raw k-mer (and its
    // reverse complement) appears in.  Using real molecule counts rather than a
    // sentinel 1 allows the graph to later filter low-evidence k-mers by
    // molecule support.
    let mut raw_kmer_counts: HashMap<u64, u32> = HashMap::new();
    for read in &reads {
        let overlaps_target = KmerIterator::new(&read.sequence, k)
            .any(|(_, km)| allowlist.contains(&canonical(km, k)));
        if overlaps_target {
            for (_, raw_km) in KmerIterator::new(&read.sequence, k) {
                *raw_kmer_counts.entry(raw_km).or_insert(0) += 1;
                *raw_kmer_counts
                    .entry(reverse_complement(raw_km, k))
                    .or_insert(0) += 1;
            }
        }
    }

    // Pass 2: build the raw graph index from the counted k-mers.
    let mut raw_graph_index = HashKmerIndex::new();
    for (raw_km, count) in raw_kmer_counts {
        raw_graph_index.insert(
            raw_km,
            kam_core::kmer::MoleculeEvidence {
                n_molecules: count,
                ..Default::default()
            },
        );
    }

    let index_qc = IndexQc {
        stage: "kmer_indexing".to_string(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        n_target_kmers,
        n_kmers_observed,
        mean_molecule_depth,
        passed: true,
    };
    write_qc(&output_dir.join("index_qc.json"), &index_qc)?;
    eprintln!(
        "[run/index] target_kmers={} observed={} time_ms={}",
        n_target_kmers,
        n_kmers_observed,
        t_index.elapsed().as_millis()
    );

    // ── Stage 3: Pathfind ─────────────────────────────────────────────────────
    let t_pathfind = std::time::Instant::now();

    // Build a single de Bruijn graph from on-target raw k-mers with sufficient
    // canonical evidence.
    //
    // raw_graph_index is a topology map (all n_molecules=1), so min_molecules
    // filtering on it is meaningless.  Instead we filter the node list against
    // the canonical evidence index, keeping only raw k-mers whose canonical
    // form is supported by ≥2 molecules.  This eliminates singleton sequencing-
    // error k-mers (which appear in exactly 1 molecule per error event) that
    // otherwise inflate the graph to 50k+ nodes and cause the DFS to explore
    // an exponentially large space of dead-end branches.
    //
    // Genuine reference k-mers survive (mean depth ~45).  Variant k-mers at
    // ≥2 supporting molecules also survive.  At very low VAF where the variant
    // is undetectable (<2 molecules), the graph simply has no alt path — the
    // correct result.
    let all_raw_kmers: Vec<u64> = raw_graph_index
        .iter()
        .filter(|(&km, _)| {
            index
                .get(canonical(km, k))
                .map(|e| e.n_molecules >= 2)
                .unwrap_or(false)
        })
        .map(|(&km, _)| km)
        .collect();
    let graph = DeBruijnGraph::from_index(&raw_graph_index, k, &all_raw_kmers, 1);

    let mut all_scored: Vec<(String, Vec<ScoredPath>)> = Vec::new();
    let mut n_targets_queried: u64 = 0;
    let mut n_targets_with_variants: u64 = 0;
    let mut n_anchors_non_unique: u64 = 0;
    let mut n_walk_budget_exceeded: u64 = 0;
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
        let alt_max_path = if cfg.input.sv_junctions.is_some() {
            target_max_path + 250
        } else {
            target_max_path
        };
        let walk_config = WalkConfig {
            max_path_length: target_max_path,
            max_expansions: 2_000_000,
            ..Default::default()
        };

        // Use evidence-biased DFS: sort successors by molecule count descending
        // so the high-evidence reference path is found before max_paths is
        // exhausted by low-evidence error-k-mer branches.
        let (raw_paths, walk_budget_hit) =
            walk_paths_biased(&graph, start_raw, end_raw, &walk_config, |raw_km| {
                index
                    .get(canonical(raw_km, k))
                    .map(|e| e.n_molecules)
                    .unwrap_or(0)
            });
        if walk_budget_hit {
            n_walk_budget_exceeded += 1;
        }

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
                // Stop the guided search once 20 distinct alt paths have been
                // found. Without this limit, the loop runs sub-DFS from every
                // high-evidence branch across all ~270 reference k-mers, even
                // after the SV alt path has been found — causing 4-minute
                // runtimes for large insertions where ~800 branches pass the
                // evidence filter.
                20,
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

            // Always try synthetic DUP and InvDel alt paths from the junction file.
            // The graph-based traversal cannot find these paths at low VAF because:
            // - DUP: the full junction k-mer chain is not in the raw graph.
            // - InvDel: the alt path starts with an inverted sequence that shares no
            //   k-mers with the reference, so the reference-guided walker cannot
            //   enter the alt branch at all.
            // Graph-derived paths (SNVs, chimeric bridges) may already be in `paths`
            // but are unrelated to the SV. Synthetic paths are added unconditionally
            // so they are scored alongside any graph-derived paths.
            if let Some(ref jdata) = sv_junctions_data {
                for (_junc_id, junc_seq) in jdata {
                    if let Some(alt_path) = synthesize_dup_alt_path(target_seq, junc_seq, k) {
                        if existing_seqs.insert(alt_path.sequence.clone()) {
                            paths.push(alt_path);
                        }
                    }
                    if let Some(alt_path) = synthesize_invdel_alt_path(target_seq, junc_seq, k) {
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
        n_targets_walk_budget_exceeded: n_walk_budget_exceeded,
        passed: true,
    };
    write_qc(&output_dir.join("pathfind_qc.json"), &pathfind_qc)?;
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
    let caller_config = cfg.to_caller_config();

    let mut all_calls = Vec::new();
    let mut n_pass: u64 = 0;
    let mut n_filtered: u64 = 0;

    // Build a map from target_id → target_seq for reference lookup.
    let target_map: HashMap<&str, &[u8]> = targets
        .iter()
        .map(|(id, seq)| (id.as_str(), seq.as_slice()))
        .collect();

    // Collect per-target reference depths for the fusion VAF denominator.
    // The reference path's mean_molecules is the best estimate of sequencing
    // depth at each target locus.
    let mut target_depths: HashMap<String, f64> = HashMap::new();
    for (target_id, scored_paths) in &all_scored {
        if let Some(ref_sp) = scored_paths.iter().find(|p| p.is_reference) {
            target_depths.insert(
                target_id.clone(),
                ref_sp.aggregate_evidence.mean_molecules as f64,
            );
        }
    }

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

    // ── Fusion target walking and calling ─────────────────────────────────────
    // Fusion targets are processed after normal targets so that partner depths
    // are available. For each fusion target: walk, score, then call using the
    // partner locus depths as the VAF denominator.
    if !fusion_targets_data.is_empty() {
        eprintln!(
            "[run/call] processing {} fusion targets",
            fusion_targets_data.len()
        );
        for ft in &fusion_targets_data {
            let fusion_seq = &ft.sequence;

            // Walk the fusion target through the graph using the same soft-anchor
            // and scoring logic as normal targets.
            let (start_offset, start_raw) = match find_soft_anchor(fusion_seq, k, &graph, false, 10)
            {
                Some(v) => v,
                None => {
                    eprintln!("[run/fusion] {}: no start anchor in graph", ft.name);
                    continue;
                }
            };
            let (end_offset, end_raw) = match find_soft_anchor(fusion_seq, k, &graph, true, 10) {
                Some(v) => v,
                None => {
                    eprintln!("[run/fusion] {}: no end anchor in graph", ft.name);
                    continue;
                }
            };

            if start_raw == end_raw {
                eprintln!(
                    "[run/fusion] {}: start and end anchors are identical, skipping",
                    ft.name
                );
                continue;
            }

            let effective_len = fusion_seq.len() - start_offset - end_offset;
            let max_path = if k > 1 {
                effective_len.saturating_sub(k - 1) + 50
            } else {
                150
            };
            let walk_config = WalkConfig {
                max_path_length: max_path,
                max_expansions: 2_000_000,
                ..Default::default()
            };

            let (raw_paths, _) =
                walk_paths_biased(&graph, start_raw, end_raw, &walk_config, |raw_km| {
                    index
                        .get(canonical(raw_km, k))
                        .map(|e| e.n_molecules)
                        .unwrap_or(0)
                });

            if raw_paths.is_empty() {
                eprintln!("[run/fusion] {}: no paths found (fusion absent)", ft.name);
                continue;
            }

            // Pad paths with the soft-anchor prefix/suffix.
            let prefix = &fusion_seq[..start_offset];
            let suffix = &fusion_seq[fusion_seq.len() - end_offset..];
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

            // Score all paths. For fusion targets, the path that matches the
            // fusion sequence is the fusion allele (flagged as is_reference by
            // score_and_rank_paths because it matches the provided reference seq).
            let scored = score_and_rank_paths(paths, &index, fusion_seq, k);

            // Use the path closest to the fusion target sequence as the fusion
            // evidence. If score_and_rank_paths flagged a reference match, use it;
            // otherwise fall back to the highest-evidence path.
            let fusion_path = scored
                .iter()
                .find(|p| p.is_reference)
                .or_else(|| scored.first());

            let Some(fp) = fusion_path else {
                eprintln!("[run/fusion] {}: scored paths empty after scoring", ft.name);
                continue;
            };

            // Look up partner depths by chromosome and locus midpoint.
            // Using the midpoint disambiguates when both partners are on the
            // same chromosome (e.g. both on chr1 in a synthetic benchmark).
            let a_mid = (ft.locus_a.start + ft.locus_a.end) / 2;
            let b_mid = (ft.locus_b.start + ft.locus_b.end) / 2;
            let partner_a_depth = find_partner_depth(&target_depths, &ft.locus_a.chrom, a_mid);
            let partner_b_depth = find_partner_depth(&target_depths, &ft.locus_b.chrom, b_mid);

            let context = FusionContext {
                partner_a_depth,
                partner_b_depth,
            };

            if let Some(fusion_call) =
                call_fusion(&fp.aggregate_evidence, &context, ft, &caller_config)
            {
                eprintln!(
                    "[run/fusion] {}: VAF={:.4} molecules={}",
                    fusion_call.name, fusion_call.vaf, fusion_call.n_molecules
                );
                // Convert FusionCall to a VariantCall for unified output.
                use kam_call::caller::{CallSource, VariantCall, VariantType};
                all_calls.push(VariantCall {
                    target_id: format!(
                        "{name}__{chrom_a}:{start_a}-{end_a}__{chrom_b}:{start_b}-{end_b}__fusion",
                        name = fusion_call.name,
                        chrom_a = fusion_call.locus_a.chrom,
                        start_a = fusion_call.locus_a.start,
                        end_a = fusion_call.locus_a.end,
                        chrom_b = fusion_call.locus_b.chrom,
                        start_b = fusion_call.locus_b.start,
                        end_b = fusion_call.locus_b.end,
                    ),
                    variant_type: VariantType::Fusion,
                    ref_sequence: fusion_seq.clone(),
                    alt_sequence: fp.path.sequence.clone(),
                    vaf: fusion_call.vaf,
                    vaf_ci_low: fusion_call.vaf_ci_low,
                    vaf_ci_high: fusion_call.vaf_ci_high,
                    n_molecules_ref: 0, // no wild-type reference path for a fusion target
                    n_molecules_alt: fusion_call.n_molecules,
                    n_duplex_alt: fusion_call.n_duplex,
                    n_simplex_alt: fusion_call.n_molecules.saturating_sub(fusion_call.n_duplex),
                    // Strand-level breakdown is not available from FusionCall; zero-fill.
                    n_simplex_fwd_alt: 0,
                    n_simplex_rev_alt: 0,
                    n_duplex_ref: 0,
                    n_simplex_ref: 0,
                    mean_alt_error_prob: 0.0,
                    min_variant_specific_duplex: fusion_call.n_duplex,
                    mean_variant_specific_molecules: fusion_call.n_molecules as f32,
                    confidence: fusion_call.confidence,
                    strand_bias_p: 1.0,
                    filter: fusion_call.filter,
                    ml_prob: None,
                    call_source: CallSource::Called,
                    rescue_min_alt_molecules: None,
                    rescue_alt_duplex: None,
                    rescue_approx_vaf: None,
                    rescue_kmers_found: None,
                    rescue_kmers_total: None,
                });
            }
        }
    }

    // Apply tumour-informed filter if target_variants is configured.
    if let Some(ref vcf_path) = cfg.input.target_variants {
        let target_set = load_target_variants(vcf_path)?;
        let tol = cfg.ti_position_tolerance();
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

    // Rescue probe: for TI targets with no matching call, query the k-mer
    // index directly for alt-supporting evidence.
    if cfg.calling.ti_rescue {
        if let Some(ref vcf_path) = cfg.input.target_variants {
            use crate::rescue::{build_alt_seq, parse_target_id, probe_ti_target};
            use kam_call::caller::{CallSource, VariantType};

            let target_set = load_target_variants(vcf_path)?;

            // Mark sub-threshold calls that match TI targets.
            for call in &mut all_calls {
                if call.filter != VariantFilter::Pass {
                    if let Some(key) = kam_call::targeting::extract_variant_key(
                        &call.target_id,
                        &call.ref_sequence,
                        &call.alt_sequence,
                    ) {
                        if target_set.contains(&key) {
                            call.call_source = CallSource::SubThreshold;
                        }
                    }
                }
            }

            // Determine which TI targets already have a PASS or SubThreshold call.
            let mut matched: std::collections::HashSet<(String, i64, String, String)> =
                std::collections::HashSet::new();
            for call in &all_calls {
                if call.filter == VariantFilter::Pass
                    || call.call_source == CallSource::SubThreshold
                {
                    if let Some(key) = kam_call::targeting::extract_variant_key(
                        &call.target_id,
                        &call.ref_sequence,
                        &call.alt_sequence,
                    ) {
                        matched.insert(key);
                    }
                }
            }

            // Probe unmatched TI targets.
            let mut rescue_calls: Vec<kam_call::caller::VariantCall> = Vec::new();

            for (chrom, vcf_pos, vcf_ref, vcf_alt) in &target_set {
                let key = (chrom.clone(), *vcf_pos, vcf_ref.clone(), vcf_alt.clone());
                if matched.contains(&key) {
                    continue;
                }

                // Find the target window that covers this position.
                let found = target_map.iter().find(|(tid, _)| {
                    if let Some((tc, ts, te)) = parse_target_id(tid) {
                        tc == chrom.as_str() && ts < *vcf_pos && *vcf_pos <= te
                    } else {
                        false
                    }
                });

                let (target_id, ref_seq) = match found {
                    Some((tid, seq)) => (*tid, *seq),
                    None => continue,
                };

                let (_, target_start, _) = match parse_target_id(target_id) {
                    Some(v) => v,
                    None => continue,
                };

                let alt_seq = match build_alt_seq(
                    ref_seq,
                    target_start,
                    *vcf_pos,
                    vcf_ref.as_bytes(),
                    vcf_alt.as_bytes(),
                ) {
                    Some(s) => s,
                    None => continue,
                };

                let evidence = probe_ti_target(&index, ref_seq, &alt_seq, k);

                let (src, min_alt, alt_dup, approx_vaf, kf, kt) = match evidence {
                    Some(ref ev) => (
                        if ev.min_alt_molecules > 0 {
                            CallSource::Rescued
                        } else {
                            CallSource::NoEvidence
                        },
                        Some(ev.min_alt_molecules),
                        Some(ev.alt_duplex),
                        Some(ev.approx_vaf),
                        Some(ev.n_alt_kmers_found),
                        Some(ev.n_alt_kmers_total),
                    ),
                    None => (
                        CallSource::NoEvidence,
                        Some(0),
                        Some(0),
                        Some(0.0),
                        Some(0),
                        Some(0),
                    ),
                };

                let variant_type = {
                    let r = vcf_ref.as_bytes();
                    let a = vcf_alt.as_bytes();
                    if r.len() == 1 && a.len() == 1 {
                        VariantType::Snv
                    } else if r.len() < a.len() {
                        VariantType::Insertion
                    } else if r.len() > a.len() {
                        VariantType::Deletion
                    } else {
                        VariantType::Mnv
                    }
                };

                let mean_ref = evidence
                    .as_ref()
                    .map(|e| e.mean_ref_molecules)
                    .unwrap_or(0.0);

                rescue_calls.push(kam_call::caller::VariantCall {
                    target_id: target_id.to_string(),
                    variant_type,
                    ref_sequence: vcf_ref.as_bytes().to_vec(),
                    alt_sequence: vcf_alt.as_bytes().to_vec(),
                    vaf: approx_vaf.unwrap_or(0.0) as f64,
                    vaf_ci_low: 0.0,
                    vaf_ci_high: 0.0,
                    n_molecules_ref: mean_ref.round() as u32,
                    n_molecules_alt: min_alt.unwrap_or(0),
                    n_duplex_alt: alt_dup.unwrap_or(0),
                    n_simplex_alt: min_alt.unwrap_or(0).saturating_sub(alt_dup.unwrap_or(0)),
                    n_simplex_fwd_alt: 0,
                    n_simplex_rev_alt: 0,
                    n_duplex_ref: 0,
                    n_simplex_ref: 0,
                    mean_alt_error_prob: 0.0,
                    min_variant_specific_duplex: alt_dup.unwrap_or(0),
                    mean_variant_specific_molecules: min_alt.unwrap_or(0) as f32,
                    confidence: 0.0,
                    strand_bias_p: 1.0,
                    filter: VariantFilter::LowConfidence,
                    ml_prob: None,
                    call_source: src,
                    rescue_min_alt_molecules: min_alt,
                    rescue_alt_duplex: alt_dup,
                    rescue_approx_vaf: approx_vaf,
                    rescue_kmers_found: kf,
                    rescue_kmers_total: kt,
                });
            }

            if !rescue_calls.is_empty() {
                eprintln!(
                    "[run/rescue] {} rescue probe records added ({} with evidence, {} with no evidence)",
                    rescue_calls.len(),
                    rescue_calls.iter().filter(|c| c.call_source == CallSource::Rescued).count(),
                    rescue_calls.iter().filter(|c| c.call_source == CallSource::NoEvidence).count(),
                );
            }
            all_calls.extend(rescue_calls);
        }
    }

    // Optional ML scoring.
    let scorer_result: Option<Result<kam_call::ml::MlScorer, _>> =
        if let Some(ref name) = args.ml_model {
            Some(crate::models::resolve(name))
        } else if let Some(ref path) = args.custom_ml_model {
            let meta_path = path.with_extension("json");
            Some(kam_call::ml::MlScorer::load(path, &meta_path))
        } else {
            None
        };

    if let Some(scorer_result) = scorer_result {
        match scorer_result {
            Ok(mut scorer) => {
                for call in &mut all_calls {
                    call.ml_prob = scorer.score(call);
                }
                eprintln!(
                    "[run/call] ML scoring applied: {} calls scored",
                    all_calls.len()
                );
            }
            Err(e) => {
                eprintln!("[run/call] WARNING: failed to load ML model: {}", e);
            }
        }
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
    write_qc(&output_dir.join("call_qc.json"), &call_qc)?;
    eprintln!(
        "[run/call] variants={} pass={} filtered={} time_ms={}",
        n_variants_called,
        n_pass,
        n_filtered,
        t_call.elapsed().as_millis()
    );

    // ── Write final variant output ────────────────────────────────────────────
    let t_output = std::time::Instant::now();
    let formats = parse_output_formats(cfg.output_format())?;
    let base_path = output_dir.join("variants");

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
        output_dir.display()
    );

    Ok(())
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Build the unified [`KamConfig`] from CLI args.
///
/// When `--config` is provided, the TOML file is loaded first, then any CLI
/// overrides are merged on top.  When `--config` is absent, the config is built
/// entirely from the CLI args.
///
/// Returns an error if the config file cannot be read or required fields are
/// absent after merging.
fn build_config(args: &RunArgs) -> Result<KamConfig, Box<dyn std::error::Error>> {
    let mut cfg = if let Some(ref config_path) = args.config {
        let mut c = KamConfig::from_file(config_path)?;
        c.merge_cli_overrides(args);
        c
    } else {
        KamConfig::from_cli(args)
    };

    // Apply hard defaults for fields that have built-in fallbacks but were not
    // supplied either by the config file or the CLI.  This ensures the config
    // is always fully specified before the pipeline starts.
    if cfg.chemistry.preset.is_none() {
        cfg.chemistry.preset = Some("twist-umi-duplex".to_string());
    }
    if cfg.chemistry.min_umi_quality.is_none() {
        cfg.chemistry.min_umi_quality = Some(20);
    }
    if cfg.assembly.min_family_size.is_none() {
        cfg.assembly.min_family_size = Some(1);
    }
    if cfg.indexing.kmer_size.is_none() {
        cfg.indexing.kmer_size = Some(31);
    }
    if cfg.output.output_format.is_none() {
        cfg.output.output_format = Some("tsv".to_string());
    }
    if cfg.calling.ti_position_tolerance.is_none() {
        cfg.calling.ti_position_tolerance = Some(0);
    }

    cfg.validate()
        .map_err(|e| -> Box<dyn std::error::Error> { e.into() })?;
    Ok(cfg)
}

/// Look up the mean molecule depth for a partner locus.
///
/// Searches `target_depths` for any entry whose key contains `chrom` as a
/// substring (e.g. "chr22" in "chr22:23632500-23632550"). Returns the mean
/// depth of the closest matching target, or 0.0 if no matching target is found.
///
/// Filters by chromosome name and, when multiple targets share the same
/// chromosome (e.g. both fusion partners on chr1), picks the one whose
/// encoded coordinate range midpoint is nearest to `locus_midpoint`.
///
/// Target IDs are expected to contain a `chrom:start-end` substring
/// (e.g. `fusion_partner_A__chr1:100-300`).  When no positional suffix is
/// parseable the function falls back to the first chromosome match.
fn find_partner_depth(
    target_depths: &HashMap<String, f64>,
    chrom: &str,
    locus_midpoint: u64,
) -> f64 {
    // Parse a `start-end` range from the `chrom:start-end` substring of `id`.
    let parse_midpoint = |id: &str| -> Option<i64> {
        // Find the last occurrence of `chrom:` in the id.
        let prefix = format!("{chrom}:");
        let pos = id.rfind(&prefix)?;
        let rest = &id[pos + prefix.len()..];
        // rest is now "start-end" possibly followed by other chars.
        let range_str = rest
            .split(|c: char| !c.is_ascii_digit() && c != '-')
            .next()?;
        let mut parts = range_str.splitn(2, '-');
        let start: i64 = parts.next()?.parse().ok()?;
        let end: i64 = parts.next()?.parse().ok()?;
        Some((start + end) / 2)
    };

    let lm = locus_midpoint as i64;
    target_depths
        .iter()
        .filter(|(id, _)| id.contains(chrom))
        .min_by_key(|(id, _)| {
            parse_midpoint(id)
                .map(|m| (m - lm).abs())
                .unwrap_or(i64::MAX)
        })
        .map(|(_, &depth)| depth)
        .unwrap_or(0.0)
}

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
/// In the VCF representation, REF = `target[dup_start]` (anchor base) and
/// ALT = `target[dup_start]` + D.  The copy is therefore inserted right after
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

/// Construct a synthetic [`GraphPath`] for an inversion-deletion (InvDel) alt allele.
///
/// InvDel alt paths cannot be found by reference-guided graph traversal because the
/// alt allele starts with an inverted sequence that shares no k-mers with the reference.
/// This function bypasses the graph by reconstructing the expected alt sequence from the
/// target and junction sequences.
///
/// The InvDel junction sequence is:
/// - Left half:  last `half` bp of the alt allele at the variant boundary.
///   This is the sequence immediately before the post-variant reference continues.
///   It is the reverse-complement of the first `half` bp of the inverted region.
/// - Right half: first `half` bp of the reference IMMEDIATELY after the variant ends.
///   This is `target[inv_end..inv_end + half]` where `inv_end` is the variant end
///   position in the target.
///
/// The alt sequence reconstructed here covers the target window as follows:
/// - If the inverted region starts within the target (`inv_start < inv_end`):
///   `target[0..inv_start] + RC(target[inv_start..inv_end]) + target[inv_end..]`
/// - If the inverted region starts before the target (no `inv_start` match):
///   `RC(target[0..inv_end - half]) + junction_left + junction_right + target[inv_end + half..]`
///
/// In the second (common) case the interior of the inverted region is approximated by
/// the RC of the reference.  The junction-spanning k-mers (those that cross the
/// inversion boundary) are reconstructed exactly from the junction sequence.  These
/// are the discriminatory k-mers that distinguish the alt allele from the reference.
///
/// The returned path uses raw (forward) k-mer encodings and is scored by the
/// caller via the canonical index.
///
/// Returns `None` if the junction does not match this target (right half absent from
/// target, or the lengths are invalid).
fn synthesize_invdel_alt_path(target_seq: &[u8], junc_seq: &[u8], k: usize) -> Option<GraphPath> {
    if junc_seq.len() < 2 || !junc_seq.len().is_multiple_of(2) {
        return None;
    }
    let half = junc_seq.len() / 2;
    if half < k {
        // Junction too short to produce any unique k-mers.
        return None;
    }

    let left_half = &junc_seq[..half];
    let right_half = &junc_seq[half..];

    // The right half of the junction is the reference sequence immediately after the
    // variant ends.  Find it in the target to locate inv_end.
    let inv_end = find_subseq(target_seq, right_half)?;

    // Try to locate where the inversion starts within the target by searching for the
    // RC of the junction left half.  The left half of the junction is the last `half`
    // bp of the alt allele, which is the RC of the last `half` bp of the inverted
    // region in the reference.  Therefore RC(left_half) = last `half` bp of the
    // inverted region as it appears in the reference (i.e. in the target).
    let rc_left = rc_seq(left_half);
    let inv_start_opt = find_subseq(target_seq, &rc_left);

    let alt_seq: Vec<u8> = if let Some(inv_start) = inv_start_opt {
        if inv_start >= inv_end {
            return None;
        }
        // The inversion is entirely within the target.
        // alt = target[0..inv_start] + RC(target[inv_start..inv_end]) + target[inv_end..]
        let inverted = rc_seq(&target_seq[inv_start..inv_end]);
        let mut seq = Vec::with_capacity(target_seq.len());
        seq.extend_from_slice(&target_seq[..inv_start]);
        seq.extend_from_slice(&inverted);
        seq.extend_from_slice(&target_seq[inv_end..]);
        seq
    } else {
        // The inversion starts before the target window (common for InvDel variants
        // where the target amplicon begins inside the inverted region).
        // Use the junction as an exact anchor for the boundary k-mers and approximate
        // the interior using the RC of the reference.
        if inv_end < half {
            return None;
        }
        let interior_rc = rc_seq(&target_seq[..inv_end - half]);
        let suffix = &target_seq[inv_end + half..];
        let mut seq = Vec::with_capacity(interior_rc.len() + junc_seq.len() + suffix.len());
        seq.extend_from_slice(&interior_rc);
        seq.extend_from_slice(left_half);
        seq.extend_from_slice(right_half);
        seq.extend_from_slice(suffix);
        seq
    };

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

/// Reverse-complement a DNA byte slice.
///
/// Non-ACGT bytes are complemented as `N` → `N`.
fn rc_seq(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        })
        .collect()
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

    /// Build a minimal [`RunArgs`] for testing, given file paths.
    ///
    /// All override fields are `None` so defaults apply.  Supply `kmer_size`
    /// via the override field when a non-default k is needed.
    fn minimal_run_args(
        r1: PathBuf,
        r2: PathBuf,
        targets: PathBuf,
        output_dir: PathBuf,
    ) -> RunArgs {
        RunArgs {
            config: None,
            r1: Some(r1),
            r2: Some(r2),
            targets: Some(targets),
            output_dir: Some(output_dir),
            chemistry_override: None,
            min_umi_quality_override: Some(0), // disable quality filter in tests
            min_family_size_override: None,
            min_template_length: None,
            kmer_size_override: None,
            min_confidence: None,
            strand_bias_threshold: None,
            min_alt_molecules: None,
            min_alt_duplex: None,
            sv_min_confidence: None,
            sv_min_alt_molecules: None,
            sv_strand_bias_threshold_override: None,
            max_vaf: None,
            sv_junctions: None,
            fusion_targets: None,
            target_variants: None,
            ti_position_tolerance_override: None,
            ti_rescue: false,
            output_format_override: None,
            qc_output: None,
            log_dir: None,
            log: vec![],
            threads: None,
            ml_model: None,
            custom_ml_model: None,
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

        let mut args = minimal_run_args(r1_path, r2_path, targets_path, output_dir.clone());
        // Use a small k so the short synthetic reads are covered.
        args.kmer_size_override = Some(8);

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

        let mut args = minimal_run_args(r1_path, r2_path, targets_path, output_dir.clone());
        args.kmer_size_override = Some(4);

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

        let mut args = minimal_run_args(r1_path, r2_path, targets_path, output_dir.clone());
        args.kmer_size_override = Some(8);
        args.output_format_override = Some("tsv,vcf".to_string());

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

    /// CLI-only mode works without a config file.
    #[test]
    fn cli_only_mode_no_config_file() {
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

        let output_dir = dir.path().join("out_cli");

        // config = None: pure CLI mode.
        let mut args = minimal_run_args(r1_path, r2_path, targets_path, output_dir.clone());
        args.kmer_size_override = Some(8);

        run_pipeline(args).expect("CLI-only mode should succeed");
        assert!(output_dir.join("variants.tsv").exists());
    }

    /// Config file values are used when CLI flags are absent.
    #[test]
    fn config_file_values_used() {
        use std::io::Write as IoWrite2;

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

        let output_dir = dir.path().join("out_cfg");

        // Write a config file supplying all required fields.
        let config_path = dir.path().join("config.toml");
        let config_content = format!(
            r#"
[input]
r1 = "{r1}"
r2 = "{r2}"
targets = "{targets}"

[output]
output_dir = "{out}"

[indexing]
kmer_size = 8

[chemistry]
min_umi_quality = 0
"#,
            r1 = r1_path.display(),
            r2 = r2_path.display(),
            targets = targets_path.display(),
            out = output_dir.display(),
        );
        {
            let mut f = std::fs::File::create(&config_path).expect("create config");
            f.write_all(config_content.as_bytes())
                .expect("write config");
        }

        // No CLI path flags — everything comes from the config file.
        let args = RunArgs {
            config: Some(config_path),
            r1: None,
            r2: None,
            targets: None,
            output_dir: None,
            chemistry_override: None,
            min_umi_quality_override: None,
            min_family_size_override: None,
            min_template_length: None,
            kmer_size_override: None,
            min_confidence: None,
            strand_bias_threshold: None,
            min_alt_molecules: None,
            min_alt_duplex: None,
            sv_min_confidence: None,
            sv_min_alt_molecules: None,
            sv_strand_bias_threshold_override: None,
            max_vaf: None,
            sv_junctions: None,
            fusion_targets: None,
            target_variants: None,
            ti_position_tolerance_override: None,
            ti_rescue: false,
            output_format_override: None,
            qc_output: None,
            log_dir: None,
            log: vec![],
            threads: None,
            ml_model: None,
            custom_ml_model: None,
        };

        run_pipeline(args).expect("config file mode should succeed");
        assert!(output_dir.join("variants.tsv").exists());
    }

    /// CLI flags override config file values.
    #[test]
    fn cli_flags_override_config_file() {
        use std::io::Write as IoWrite2;

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

        // Config file points to a non-existent output dir (will be overridden).
        let config_output_dir = dir.path().join("cfg_out");
        // CLI overrides to a different output dir.
        let cli_output_dir = dir.path().join("cli_out");

        let config_path = dir.path().join("config.toml");
        let config_content = format!(
            r#"
[input]
r1 = "{r1}"
r2 = "{r2}"
targets = "{targets}"

[output]
output_dir = "{out}"
output_format = "vcf"

[indexing]
kmer_size = 4

[chemistry]
min_umi_quality = 0
"#,
            r1 = r1_path.display(),
            r2 = r2_path.display(),
            targets = targets_path.display(),
            out = config_output_dir.display(),
        );
        {
            let mut f = std::fs::File::create(&config_path).expect("create config");
            f.write_all(config_content.as_bytes())
                .expect("write config");
        }

        let args = RunArgs {
            config: Some(config_path),
            r1: None,
            r2: None,
            targets: None,
            // CLI overrides output_dir and output_format over config values.
            output_dir: Some(cli_output_dir.clone()),
            chemistry_override: None,
            min_umi_quality_override: None,
            min_family_size_override: None,
            min_template_length: None,
            kmer_size_override: Some(8),
            min_confidence: None,
            strand_bias_threshold: None,
            min_alt_molecules: None,
            min_alt_duplex: None,
            sv_min_confidence: None,
            sv_min_alt_molecules: None,
            sv_strand_bias_threshold_override: None,
            max_vaf: None,
            sv_junctions: None,
            fusion_targets: None,
            target_variants: None,
            ti_position_tolerance_override: None,
            ti_rescue: false,
            output_format_override: Some("tsv".to_string()),
            qc_output: None,
            log_dir: None,
            log: vec![],
            threads: None,
            ml_model: None,
            custom_ml_model: None,
        };

        run_pipeline(args).expect("CLI override mode should succeed");

        // Output should be in the CLI-specified dir, not the config dir.
        assert!(
            cli_output_dir.join("variants.tsv").exists(),
            "output should be in CLI-specified dir"
        );
        assert!(
            !config_output_dir.exists(),
            "config-specified dir should not be created"
        );
    }

    // ── synthesize_invdel_alt_path unit tests ─────────────────────────────────

    /// Inversion fully within the target: inv_start and inv_end both found.
    ///
    /// Layout:
    ///   target = AAAACCCCTTTT (12bp)
    ///   inverted region = target[4..8] = CCCC
    ///   inv_start = 4 (RC(junction_left) = RC(GGGG) = CCCC found at target[4])
    ///   inv_end   = 8 (junction_right = TTTT found at target[8])
    ///   alt = AAAA + RC(CCCC) + TTTT = AAAAGGGGT TTT → AAAAGGGGTTTT
    #[test]
    fn synthesize_invdel_inversion_fully_within_target() {
        // target:       AAAA CCCC TTTT
        //                    ^^^^ inverted region (target[4..8])
        // alt:          AAAA GGGG TTTT
        //
        // junction_left = last 4bp of alt at the boundary = GGGG
        //   (= RC of CCCC = last half of the inverted region)
        // RC(junction_left) = RC(GGGG) = CCCC → found at target[4] = inv_start
        // junction_right = TTTT → found at target[8] = inv_end
        let target = b"AAAACCCCTTTT"; // 12bp
        let junc = b"GGGGTTTT"; // half=4; left=GGGG, right=TTTT
        let k = 4;

        let path = synthesize_invdel_alt_path(target, junc, k).expect("should produce a path");

        let expected: &[u8] = b"AAAAGGGGTTTT";
        assert_eq!(
            path.sequence, expected,
            "fully-within inversion should be AAAA+GGGG+TTTT"
        );
        assert!(!path.kmers.is_empty(), "path should have k-mers");
        assert_eq!(path.length, path.kmers.len(), "length field must match");
    }

    /// Inversion starting before the target window (the common InvDel case).
    ///
    /// The RC of junction_left does not appear in the target (because the start of
    /// the inverted region is before the target).  The fallback path is:
    ///   RC(target[0..inv_end - half]) + junction_left + junction_right + target[inv_end + half..]
    ///
    /// Layout:
    ///   target = TTTTACGTACGT (12bp)
    ///   inv_end = 4 (ACGT starts at target[4])
    ///   half = 4
    ///   junction_left = CCCC (RC=GGGG, not in target → inversion started before target)
    ///   interior_rc = RC(target[0..0]) = "" (empty, since inv_end - half = 0)
    ///   alt = "" + CCCC + ACGT + target[8..] = CCCC + ACGT + ACGT = CCCCACGTACGT
    #[test]
    fn synthesize_invdel_inversion_starts_before_target() {
        let target = b"TTTTACGTACGT";
        // junction: left=CCCC (RC=GGGG, absent from target → uses fallback)
        //           right=ACGT (appears at target[4] → inv_end=4)
        let junc = b"CCCCACGT"; // 8 bytes, half=4
        let k = 4;

        let path = synthesize_invdel_alt_path(target, junc, k)
            .expect("fallback path should be produced for inversion-before-target");

        // inv_end=4, half=4, interior_rc = RC(target[0..0]) = ""
        // alt = "" + CCCC + ACGT + target[8..12] = CCCCACGTACGT
        let expected: Vec<u8> = b"CCCCACGTACGT".to_vec();
        assert_eq!(path.sequence, expected);
        assert!(!path.kmers.is_empty());
    }

    /// inv_end < half causes None because the interior slice would underflow.
    #[test]
    fn synthesize_invdel_inv_end_less_than_half_returns_none() {
        // target = ACGTTTTT (8bp)
        // right_half = ACGT appears at position 0. inv_end=0, half=4. 0 < 4 → None.
        let target = b"ACGTTTTT";
        let junc = b"CCCCACGT"; // right half ACGT found at pos 0, inv_end=0 < half=4
        assert!(
            synthesize_invdel_alt_path(target, junc, 4).is_none(),
            "inv_end < half should return None to prevent underflow"
        );
    }

    /// Junction with odd length returns None.
    #[test]
    fn synthesize_invdel_odd_junction_returns_none() {
        let target = b"ACGTACGT";
        let junc = b"ACGTACG"; // 7 bytes — odd, not multiple of 2
        assert!(
            synthesize_invdel_alt_path(target, junc, 4).is_none(),
            "odd-length junction should return None"
        );
    }

    /// Right half not found in target returns None.
    #[test]
    fn synthesize_invdel_right_half_absent_returns_none() {
        let target = b"ACGTACGT";
        // right half = TTTT, not in target
        let junc = b"CCCCTTTT";
        assert!(
            synthesize_invdel_alt_path(target, junc, 4).is_none(),
            "absent right half should return None"
        );
    }

    /// Junction shorter than 2*k returns None (no discriminatory k-mers possible).
    #[test]
    fn synthesize_invdel_junction_too_short_returns_none() {
        let target = b"ACGTACGT";
        // junction is 4 bytes (half=2), k=4, half < k → None
        let junc = b"ACGT";
        assert!(
            synthesize_invdel_alt_path(target, junc, 4).is_none(),
            "junction shorter than 2*k should return None"
        );
    }

    /// rc_seq correctly reverse-complements a simple sequence.
    #[test]
    fn rc_seq_correctness() {
        assert_eq!(rc_seq(b"ACGT"), b"ACGT"); // ACGT is its own RC
        assert_eq!(rc_seq(b"AAAA"), b"TTTT");
        assert_eq!(rc_seq(b"CCCC"), b"GGGG");
        assert_eq!(rc_seq(b"ATCG"), b"CGAT");
        assert_eq!(rc_seq(b""), b"" as &[u8]);
    }
}
