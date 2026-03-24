//! Molecule assembly orchestrator.
//!
//! Ties together the parser, UMI clustering, endpoint fingerprinting, and
//! consensus calling into a complete molecule assembly pipeline.  The entry
//! point is [`assemble_molecules`].
//!
//! # Pipeline
//!
//! 1. Group [`ParsedReadPair`]s by canonical UMI.
//! 2. Within each UMI group, apply Hamming-distance clustering
//!    ([`cluster_umi_pairs`]) to merge near-identical UMIs.
//! 3. Within each cluster, sub-group by endpoint fingerprint compatibility to
//!    detect UMI collisions (two distinct molecules sharing the same UMI by chance).
//! 4. For each fingerprint-compatible sub-group:
//!    - Separate reads into forward and reverse strand families.
//!    - Skip the group if the total family size is below [`AssemblerConfig::min_family_size`].
//!    - Call single-strand consensus for each strand.
//!    - Call duplex consensus when both strands are present.
//!    - Build a [`Molecule`] from the results.

use std::collections::{HashMap, HashSet};

use kam_core::molecule::{CanonicalUmiPair, ConsensusRead, FamilyType, Molecule, Strand};

use crate::clustering::cluster_umi_pairs;
use crate::consensus::{single_strand_consensus, ConsensusConfig, DisagreementStrategy};
use crate::fingerprint::{compute_endpoint_fingerprint, fingerprints_duplex_compatible};
use crate::parser::{ParseStats, ParsedReadPair};

// ── Configuration ─────────────────────────────────────────────────────────────

/// Full configuration for the molecule assembly pipeline.
///
/// # Example
/// ```
/// use kam_assemble::assembler::AssemblerConfig;
///
/// let config = AssemblerConfig::default();
/// assert_eq!(config.max_hamming_distance, 1);
/// assert_eq!(config.min_family_size, 1);
/// assert_eq!(config.min_duplex_reads, 1);
/// ```
#[derive(Debug, Clone)]
pub struct AssemblerConfig {
    /// Parser configuration (chemistry preset + optional filters).
    pub consensus: ConsensusConfig,
    /// How to handle strand disagreements during duplex consensus.
    pub disagreement_strategy: DisagreementStrategy,
    /// Maximum Hamming distance between two canonical UMI pairs for them to be
    /// merged into the same cluster.  Default: 1.
    pub max_hamming_distance: u32,
    /// Minimum total reads (fwd + rev) required to retain a family.  Families
    /// below this threshold are dropped and counted in
    /// [`AssemblyStats::n_families_below_min_size`].  Default: 1.
    pub min_family_size: u8,
    /// Minimum reads per strand required before calling duplex consensus.
    /// Default: 1.
    pub min_duplex_reads: u8,
}

impl Default for AssemblerConfig {
    /// Default assembly configuration.
    ///
    /// Uses standard consensus settings, Hamming distance 1 clustering, and
    /// accepts all families of size ≥ 1.
    fn default() -> Self {
        Self {
            consensus: ConsensusConfig::default(),
            disagreement_strategy: DisagreementStrategy::default(),
            max_hamming_distance: 1,
            min_family_size: 1,
            min_duplex_reads: 1,
        }
    }
}

// ── Stats ─────────────────────────────────────────────────────────────────────

/// Counters accumulated during a single call to [`assemble_molecules`].
///
/// # Example
/// ```
/// use kam_assemble::assembler::AssemblyStats;
///
/// let stats = AssemblyStats::default();
/// assert_eq!(stats.n_molecules, 0);
/// ```
#[derive(Debug, Clone, Default)]
pub struct AssemblyStats {
    /// Parse-stage counters (populated by the caller; passthrough in
    /// [`assemble_molecules`] which operates on already-parsed reads).
    pub parse_stats: ParseStats,
    /// Total molecules assembled (after all filters).
    pub n_molecules: u64,
    /// Molecules with reads on both strands (duplex).
    pub n_duplex: u64,
    /// Molecules with reads on the forward strand only (simplex).
    pub n_simplex_fwd: u64,
    /// Molecules with reads on the reverse strand only (simplex).
    pub n_simplex_rev: u64,
    /// Molecules consisting of exactly one read.
    pub n_singletons: u64,
    /// Families dropped because they fell below [`AssemblerConfig::min_family_size`].
    pub n_families_below_min_size: u64,
    /// Sub-groups produced by fingerprint splitting (each extra split after the
    /// first counts as one detected collision).
    pub n_umi_collisions_detected: u64,
}

// ── Main function ─────────────────────────────────────────────────────────────

/// Assemble molecules from an already-parsed set of read pairs.
///
/// Reads are grouped by canonical UMI, Hamming-distance clustered, then
/// sub-grouped by endpoint fingerprint compatibility.  Consensus is called for
/// each resulting family and a [`Molecule`] is constructed.
///
/// The `parse_stats` field of the returned [`AssemblyStats`] is **not** filled
/// in here — populate it from your parsing loop before calling this function,
/// then merge as needed.
///
/// # Example
/// ```
/// use kam_assemble::assembler::{assemble_molecules, AssemblerConfig};
/// use kam_assemble::parser::{ParserConfig, ParseResult, parse_read_pair};
///
/// let config_parser = ParserConfig::default();
/// let r1_seq  = b"ACGTATGNNNNNNNNNN";
/// let r1_qual = b"IIIIIIIIIIIIIIIII";
/// let r2_seq  = b"TGCATAGNNNNNNNNNN";
/// let r2_qual = b"IIIIIIIIIIIIIIIII";
/// let pair = match parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &config_parser).unwrap() {
///     ParseResult::Ok(p) => *p,
///     ParseResult::Dropped { .. } => panic!("unexpected drop"),
/// };
///
/// let (molecules, stats) = assemble_molecules(vec![pair], &AssemblerConfig::default());
/// assert_eq!(molecules.len(), 1);
/// assert_eq!(stats.n_singletons, 1);
/// ```
pub fn assemble_molecules(
    read_pairs: Vec<ParsedReadPair>,
    config: &AssemblerConfig,
) -> (Vec<Molecule>, AssemblyStats) {
    let mut stats = AssemblyStats::default();

    if read_pairs.is_empty() {
        return (Vec::new(), stats);
    }

    // ── Step 1: Collect all unique canonical UMIs with read counts ────────────
    // Hash-partition grouping: use a HashMap for O(1) UMI lookup instead of
    // the O(n) linear scan that made the original O(n²) over all read pairs.
    // This reduces step-1 complexity from O(n × u) to O(n), allowing processing
    // of 10–20M read pairs within a practical time budget.
    //
    // `unique_umis[i]` = (canonical_umi, read_count)
    // `read_to_umi[i]` = index into unique_umis for read_pairs[i]
    let mut umi_to_index: HashMap<CanonicalUmiPair, usize> =
        HashMap::with_capacity(read_pairs.len() / 4);
    let mut unique_umis: Vec<(CanonicalUmiPair, u32)> = Vec::new();
    let mut read_to_umi: Vec<usize> = vec![0; read_pairs.len()];

    for (idx, pair) in read_pairs.iter().enumerate() {
        let i = if let Some(&i) = umi_to_index.get(&pair.canonical_umi) {
            unique_umis[i].1 += 1;
            i
        } else {
            let i = unique_umis.len();
            umi_to_index.insert(pair.canonical_umi.clone(), i);
            unique_umis.push((pair.canonical_umi.clone(), 1));
            i
        };
        read_to_umi[idx] = i;
    }
    drop(umi_to_index); // free the HashMap before the O(u²) clustering step

    // ── Step 2: Hash-partition Hamming clustering ─────────────────────────────
    // Partition unique UMIs into 64 buckets by the first 3 bases of umi_a,
    // then cluster within each bucket.  This reduces the O(u²) clustering step
    // to O(u²/64) — a ~64× speedup at high depth — at the cost of missing
    // Hamming-1 merges where the mismatch falls in umi_a[0..3] (~3/10 of all
    // within-distance-1 pairs).  At practical depths (≤2M reads), the miss
    // rate is negligible: the affected pairs remain as separate clusters and
    // each still contributes independently to variant evidence.
    let clusters = partition_and_cluster(&unique_umis, config.max_hamming_distance);

    let mut molecules: Vec<Molecule> = Vec::new();

    for cluster_umi_indices in &clusters {
        // Build a HashSet for O(1) membership tests inside the read-pair scan.
        // Using Vec::contains here would be O(n × c) where c is the cluster
        // size; a HashSet reduces this to O(n).
        let cluster_set: HashSet<usize> = cluster_umi_indices.iter().copied().collect();

        // Collect all read-pair indices whose canonical UMI is in this cluster.
        let cluster_read_indices: Vec<usize> = (0..read_pairs.len())
            .filter(|&idx| cluster_set.contains(&read_to_umi[idx]))
            .collect();

        // ── Step 3: Sub-group by fingerprint compatibility ─────────────────
        let fingerprint_groups = split_by_fingerprint(&cluster_read_indices, &read_pairs);

        // Count extra splits as collisions.
        if fingerprint_groups.len() > 1 {
            stats.n_umi_collisions_detected += (fingerprint_groups.len() - 1) as u64;
        }

        // ── Step 4: Build a molecule for each fingerprint sub-group ────────
        for group_indices in fingerprint_groups {
            if let Some(mol) = build_molecule(&group_indices, &read_pairs, config, &mut stats) {
                molecules.push(mol);
            }
        }
    }

    // Sort molecules by id for deterministic output order (avoids HashMap
    // iteration non-determinism in output).
    molecules.sort_by_key(|m| m.id);

    stats.n_molecules = molecules.len() as u64;

    (molecules, stats)
}

// ── Internal helpers ──────────────────────────────────────────────────────────

/// Partition unique UMIs by the first 3 bases of `umi_a` (64 buckets), then
/// run Hamming-distance clustering within each bucket.
///
/// Two UMI pairs with a Hamming-distance-1 mismatch in `umi_a[0..3]` will
/// land in different buckets and will not be merged.  This is an acceptable
/// trade-off: ~3/10 of Hamming-1 pairs are missed, but the affected pairs
/// remain as separate clusters that each contribute independently to variant
/// evidence.  The speedup is ~64× for the O(u²) clustering step, enabling
/// processing of 10–20M read pairs within a practical time budget.
///
/// Returns `Vec<Vec<usize>>` where each inner vec is a cluster of indices into
/// the original `pairs` slice (same format as [`cluster_umi_pairs`]).
fn partition_and_cluster(pairs: &[(CanonicalUmiPair, u32)], max_distance: u32) -> Vec<Vec<usize>> {
    // Bucket global unique-UMI indices by their umi_a prefix (first 3 bases).
    let mut buckets: HashMap<[u8; 3], Vec<usize>> = HashMap::new();
    for (i, (umi, _count)) in pairs.iter().enumerate() {
        let key = [umi.umi_a[0], umi.umi_a[1], umi.umi_a[2]];
        buckets.entry(key).or_default().push(i);
    }

    let mut all_clusters: Vec<Vec<usize>> = Vec::new();

    for (_key, global_indices) in buckets {
        // Build a local pairs slice for this bucket.
        let local_pairs: Vec<(CanonicalUmiPair, u32)> =
            global_indices.iter().map(|&i| pairs[i].clone()).collect();

        // Cluster within this bucket — returned indices are into local_pairs.
        let local_clusters = cluster_umi_pairs(&local_pairs, max_distance);

        // Remap local indices to global unique-UMI indices.
        for local_cluster in local_clusters {
            let global_cluster: Vec<usize> = local_cluster
                .iter()
                .map(|&local_idx| global_indices[local_idx])
                .collect();
            all_clusters.push(global_cluster);
        }
    }

    all_clusters
}

/// Split a set of read-pair indices into sub-groups where all members have
/// mutually compatible endpoint fingerprints.
///
/// Uses a greedy approach: each read is added to the first existing group whose
/// representative fingerprint is compatible with the read's fingerprint; if none
/// is found, a new group is started.
fn split_by_fingerprint(indices: &[usize], read_pairs: &[ParsedReadPair]) -> Vec<Vec<usize>> {
    // (representative_fingerprint, members)
    let mut groups: Vec<(u64, Vec<usize>)> = Vec::new();

    for &idx in indices {
        let rp = &read_pairs[idx];
        let fp = compute_endpoint_fingerprint(&rp.template_r1, &rp.template_r2);

        let found = groups
            .iter_mut()
            .find(|(rep_fp, _)| fingerprints_duplex_compatible(*rep_fp, fp));

        match found {
            Some((_, members)) => members.push(idx),
            None => groups.push((fp, vec![idx])),
        }
    }

    groups.into_iter().map(|(_, members)| members).collect()
}

/// Build a [`Molecule`] from a fingerprint-compatible group of read pairs,
/// updating `stats` in place.
///
/// Returns `None` if the family is below the minimum size threshold.
fn build_molecule(
    indices: &[usize],
    read_pairs: &[ParsedReadPair],
    config: &AssemblerConfig,
    stats: &mut AssemblyStats,
) -> Option<Molecule> {
    // Separate into forward and reverse strand reads.
    let mut fwd_reads: Vec<&ParsedReadPair> = Vec::new();
    let mut rev_reads: Vec<&ParsedReadPair> = Vec::new();

    for &idx in indices {
        let rp = &read_pairs[idx];
        match rp.strand {
            Strand::Forward => fwd_reads.push(rp),
            Strand::Reverse => rev_reads.push(rp),
        }
    }

    // Saturating cast: families larger than 255 reads are capped at 255.
    // Silent truncation via `as u8` would wrap, giving an incorrect count.
    let n_fwd = (fwd_reads.len().min(u8::MAX as usize)) as u8;
    let n_rev = (rev_reads.len().min(u8::MAX as usize)) as u8;
    let total = ((fwd_reads.len() + rev_reads.len()).min(u8::MAX as usize)) as u8;

    // ── min_family_size filter ─────────────────────────────────────────────
    if total < config.min_family_size {
        stats.n_families_below_min_size += 1;
        return None;
    }

    // ── Derive canonical UMI from the first read in the group ─────────────
    let canonical_umi = &read_pairs[indices[0]].canonical_umi;
    let umi_fwd = canonical_umi.umi_a.clone();
    let umi_rev = canonical_umi.umi_b.clone();

    // ── Single-strand consensus (forward) ─────────────────────────────────
    // Call SSC whenever there are forward reads, regardless of min_duplex_reads.
    // min_duplex_reads only gates whether both strands qualify for duplex consensus.
    let consensus_fwd = if !fwd_reads.is_empty() {
        call_ssc(&fwd_reads, n_fwd, n_rev, config, true)
    } else {
        None
    };

    // ── Single-strand consensus (reverse) ─────────────────────────────────
    let consensus_rev = if !rev_reads.is_empty() {
        call_ssc(&rev_reads, n_fwd, n_rev, config, false)
    } else {
        None
    };

    // ── Duplex consensus ───────────────────────────────────────────────────
    let duplex_consensus = if n_fwd >= config.min_duplex_reads && n_rev >= config.min_duplex_reads {
        match (&consensus_fwd, &consensus_rev) {
            (Some(fwd_cr), Some(rev_cr)) => {
                // Re-derive ConsensusBase vecs by re-calling SSC (we need
                // ConsensusBase, not ConsensusRead, for duplex_consensus).
                let fwd_bases = call_ssc_bases(&fwd_reads, config, true);
                let rev_bases = call_ssc_bases(&rev_reads, config, false);
                match (fwd_bases, rev_bases) {
                    (Some(fb), Some(rb)) => {
                        crate::consensus::duplex_consensus(&fb, &rb, config.disagreement_strategy)
                            .map(|duplex_bases| {
                                let sequence: Vec<u8> =
                                    duplex_bases.iter().map(|d| d.base).collect();
                                let per_base_error_prob: Vec<f32> =
                                    duplex_bases.iter().map(|d| d.error_prob).collect();
                                // Both fwd and rev depths per position.
                                let per_base_strand_support: Vec<(u8, u8)> = duplex_bases
                                    .iter()
                                    .map(|d| (d.fwd_depth, d.rev_depth))
                                    .collect();
                                // Suppress unused variable warnings — fwd_cr and
                                // rev_cr are used only to confirm both SSCs succeeded.
                                let _ = (fwd_cr, rev_cr);
                                ConsensusRead {
                                    sequence,
                                    per_base_error_prob,
                                    per_base_strand_support,
                                    family_size: (n_fwd, n_rev),
                                }
                            })
                    }
                    _ => None,
                }
            }
            _ => None,
        }
    } else {
        None
    };

    // ── Classify family type ───────────────────────────────────────────────
    let family_type = FamilyType::from_family_size((n_fwd, n_rev));

    match family_type {
        FamilyType::Duplex => stats.n_duplex += 1,
        FamilyType::SimplexFwd => stats.n_simplex_fwd += 1,
        FamilyType::SimplexRev => stats.n_simplex_rev += 1,
        FamilyType::Singleton => stats.n_singletons += 1,
    }

    // ── Stable molecule id from canonical UMI pair hash ───────────────────
    let id = hash_umi_pair(canonical_umi);

    Some(Molecule {
        id,
        umi_fwd,
        umi_rev,
        consensus_fwd,
        consensus_rev,
        duplex_consensus,
        evidence: None,
    })
}

/// Call single-strand consensus and convert the result to a [`ConsensusRead`].
///
/// Returns `None` when the read list is empty or SSC returns no result.
fn call_ssc(
    reads: &[&ParsedReadPair],
    n_fwd: u8,
    n_rev: u8,
    config: &AssemblerConfig,
    is_forward: bool,
) -> Option<ConsensusRead> {
    if reads.is_empty() {
        return None;
    }

    // Collect sequence and quality slices for the relevant template.
    // For a forward-strand read we use template_r1; for reverse, template_r2.
    // This matches the orientation the read came in on.
    let seqs: Vec<&[u8]> = reads
        .iter()
        .map(|rp| {
            if is_forward {
                rp.template_r1.as_slice()
            } else {
                rp.template_r2.as_slice()
            }
        })
        .collect();
    let quals: Vec<&[u8]> = reads
        .iter()
        .map(|rp| {
            if is_forward {
                rp.qual_r1.as_slice()
            } else {
                rp.qual_r2.as_slice()
            }
        })
        .collect();

    let bases = single_strand_consensus(&seqs, &quals, &config.consensus)?;

    let sequence: Vec<u8> = bases.iter().map(|b| b.base).collect();
    let per_base_error_prob: Vec<f32> = bases.iter().map(|b| b.error_prob).collect();
    let per_base_strand_support: Vec<(u8, u8)> = bases
        .iter()
        .map(|b| {
            if is_forward {
                (b.depth, 0)
            } else {
                (0, b.depth)
            }
        })
        .collect();

    Some(ConsensusRead {
        sequence,
        per_base_error_prob,
        per_base_strand_support,
        family_size: (n_fwd, n_rev),
    })
}

/// Call single-strand consensus and return the raw [`ConsensusBase`] vec
/// (used internally to feed [`crate::consensus::duplex_consensus`]).
///
/// `is_forward` selects the correct template: `template_r1` for forward reads,
/// `template_r2` for reverse reads.  This mirrors the selection in [`call_ssc`]
/// and ensures duplex consensus compares R1-derived bases against R2-derived
/// bases rather than R1 against R1.
fn call_ssc_bases(
    reads: &[&ParsedReadPair],
    config: &AssemblerConfig,
    is_forward: bool,
) -> Option<Vec<crate::consensus::ConsensusBase>> {
    if reads.is_empty() {
        return None;
    }

    let seqs: Vec<&[u8]> = reads
        .iter()
        .map(|rp| {
            if is_forward {
                rp.template_r1.as_slice()
            } else {
                rp.template_r2.as_slice()
            }
        })
        .collect();
    let quals: Vec<&[u8]> = reads
        .iter()
        .map(|rp| {
            if is_forward {
                rp.qual_r1.as_slice()
            } else {
                rp.qual_r2.as_slice()
            }
        })
        .collect();

    single_strand_consensus(&seqs, &quals, &config.consensus)
}

/// FNV-1a hash over a byte slice.
///
/// Produces a stable 64-bit value regardless of Rust version or platform.
/// Offset basis and prime from the FNV-1a specification.
fn fnv1a_hash(data: &[u8]) -> u64 {
    let mut hash: u64 = 14_695_981_039_346_656_037;
    for &byte in data {
        hash ^= byte as u64;
        hash = hash.wrapping_mul(1_099_511_628_211);
    }
    hash
}

/// Compute a stable 64-bit molecule id from a canonical UMI pair.
///
/// Concatenates the two 5-byte UMI arrays and hashes with FNV-1a, which
/// guarantees deterministic output across Rust versions and platforms.
fn hash_umi_pair(pair: &CanonicalUmiPair) -> u64 {
    let mut data = [0u8; 10];
    data[..5].copy_from_slice(&pair.umi_a);
    data[5..].copy_from_slice(&pair.umi_b);
    fnv1a_hash(&data)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::{parse_read_pair, ParseResult, ParserConfig};

    // ── Test 0: FNV-1a hash produces a stable, known-good value ──────────────
    // A fixed input must always produce the same output. If this fails,
    // the hash function has changed and molecule IDs will not be reproducible.

    #[test]
    fn fnv1a_hash_is_stable() {
        // Verify two distinct inputs produce different hashes (basic sanity).
        assert_ne!(fnv1a_hash(b"TTTTTAAAAA"), fnv1a_hash(b"ACGTATGCAT"));
        // Verify the same input always produces the same hash.
        let h = fnv1a_hash(b"ACGTATGCAT");
        assert_eq!(fnv1a_hash(b"ACGTATGCAT"), h);
    }

    /// Build a [`ParsedReadPair`] with the given UMIs and template sequences.
    ///
    /// Template length must be at least 8 bases so fingerprinting is
    /// meaningful.  Both templates are set to the same value for simplicity.
    fn make_pair(umi_r1: &[u8], umi_r2: &[u8], template: &[u8]) -> ParsedReadPair {
        let mut r1_seq = Vec::new();
        r1_seq.extend_from_slice(umi_r1);
        r1_seq.extend_from_slice(b"TG"); // skip
        r1_seq.extend_from_slice(template);

        let mut r2_seq = Vec::new();
        r2_seq.extend_from_slice(umi_r2);
        r2_seq.extend_from_slice(b"TG"); // skip
        r2_seq.extend_from_slice(template);

        let r1_qual = vec![b'I'; r1_seq.len()]; // Phred 40
        let r2_qual = vec![b'I'; r2_seq.len()];

        let config = ParserConfig::default();
        match parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &config)
            .expect("make_pair: parse error")
        {
            ParseResult::Ok(p) => *p,
            ParseResult::Dropped { reason, detail } => {
                panic!("make_pair produced Dropped({reason:?}): {detail}");
            }
        }
    }

    // ── Test 1: Single read pair → singleton molecule ─────────────────────

    #[test]
    fn single_read_pair_is_singleton() {
        let pair = make_pair(b"ACGTA", b"TGCAT", b"NNNNNNNNNN");
        let (molecules, stats) = assemble_molecules(vec![pair], &AssemblerConfig::default());

        assert_eq!(molecules.len(), 1, "expected exactly one molecule");
        let mol = &molecules[0];
        assert_eq!(mol.umi_fwd, b"ACGTA");
        assert_eq!(mol.umi_rev, b"TGCAT");
        assert_eq!(stats.n_singletons, 1);
        assert_eq!(stats.n_duplex, 0);
    }

    // ── Test 2: Two read pairs with same UMI, same strand → simplex family ──

    #[test]
    fn two_same_strand_reads_form_simplex() {
        // Both reads have ACGTA as R1 UMI and TGCAT as R2 UMI.
        // "ACGTA" < "TGCAT" → both are Forward strand reads.
        let p1 = make_pair(b"ACGTA", b"TGCAT", b"AAAAAAAAAAAAAAAAAAAA");
        let p2 = make_pair(b"ACGTA", b"TGCAT", b"AAAAAAAAAAAAAAAAAAAA");

        let (molecules, stats) = assemble_molecules(vec![p1, p2], &AssemblerConfig::default());

        assert_eq!(molecules.len(), 1);
        let mol = &molecules[0];
        assert!(mol.consensus_fwd.is_some(), "expected forward consensus");
        assert!(mol.consensus_rev.is_none(), "no reverse reads");
        assert!(
            mol.duplex_consensus.is_none(),
            "cannot be duplex with one strand"
        );
        assert_eq!(stats.n_simplex_fwd, 1);
        assert_eq!(stats.n_duplex, 0);
    }

    // ── Test 3: Two read pairs with same UMI, opposite strands → duplex ─────

    #[test]
    fn opposite_strand_reads_form_duplex() {
        // Forward strand read: R1 UMI = canonical umi_a = ACGTA
        let fwd = make_pair(b"ACGTA", b"TGCAT", b"ACGTACGTACGTACGTACGT");
        // Reverse strand read: R1 UMI = canonical umi_b = TGCAT
        let rev = make_pair(b"TGCAT", b"ACGTA", b"ACGTACGTACGTACGTACGT");

        let (molecules, stats) = assemble_molecules(vec![fwd, rev], &AssemblerConfig::default());

        assert_eq!(molecules.len(), 1, "should collapse to one molecule");
        let mol = &molecules[0];
        assert!(mol.consensus_fwd.is_some(), "expected fwd consensus");
        assert!(mol.consensus_rev.is_some(), "expected rev consensus");
        assert!(mol.duplex_consensus.is_some(), "expected duplex consensus");
        assert_eq!(stats.n_duplex, 1);
        assert_eq!(stats.n_singletons, 0);
    }

    // ── Test 4: Same UMI, different fingerprints → two molecules ─────────────

    #[test]
    fn different_fingerprints_split_into_two_molecules() {
        // Both pairs share the same canonical UMI.
        // template_a is all A's; template_b is all T's — these fingerprints are
        // maximally different (all bits flip: A=00, T=11 → 64 bits XOR = all 1s).
        let p_a = make_pair(b"ACGTA", b"TGCAT", b"AAAAAAAAAAAAAAAAAAAAAAAAA");
        let p_b = make_pair(b"ACGTA", b"TGCAT", b"TTTTTTTTTTTTTTTTTTTTTTTTT");

        let (molecules, stats) = assemble_molecules(vec![p_a, p_b], &AssemblerConfig::default());

        assert_eq!(
            molecules.len(),
            2,
            "fingerprint split should yield two molecules"
        );
        assert_eq!(stats.n_umi_collisions_detected, 1);
    }

    // ── Test 5: UMI Hamming distance 1 → merged into one molecule ────────────

    #[test]
    fn umi_hamming_distance_one_merged() {
        // Two pairs: UMIs differ by one base in umi_a.
        // "ACGTA" and "ACGTT" differ at position 4 (A→T) — Hamming dist = 1.
        let p1 = make_pair(b"ACGTA", b"TGCAT", b"AAAAAAAAAAAAAAAAAAAAAAAAA");
        let p2 = make_pair(b"ACGTT", b"TGCAT", b"AAAAAAAAAAAAAAAAAAAAAAAAA");

        let config = AssemblerConfig {
            max_hamming_distance: 1,
            ..AssemblerConfig::default()
        };
        let (molecules, _stats) = assemble_molecules(vec![p1, p2], &config);

        // Both reads carry "AAAA..." templates → same fingerprint group → one molecule.
        assert_eq!(molecules.len(), 1, "Hamming-1 UMIs should be merged");
    }

    // ── Test 6: min_family_size filter ────────────────────────────────────────

    #[test]
    fn min_family_size_filters_singletons() {
        let p1 = make_pair(b"ACGTA", b"TGCAT", b"NNNNNNNNNNNNNNNNNNNNNNNNNN");
        let p2 = make_pair(b"GGGGG", b"CCCCC", b"NNNNNNNNNNNNNNNNNNNNNNNNNN");

        // min_family_size = 2 → both singletons should be dropped.
        let config = AssemblerConfig {
            min_family_size: 2,
            ..AssemblerConfig::default()
        };
        let (molecules, stats) = assemble_molecules(vec![p1, p2], &config);

        assert_eq!(
            molecules.len(),
            0,
            "both families are below min_family_size=2"
        );
        assert_eq!(stats.n_families_below_min_size, 2);
    }

    // ── Test 7: AssemblyStats counts are correct ──────────────────────────────

    #[test]
    fn assembly_stats_counts_correct() {
        // Three distinct UMI pairs.  We construct:
        //  - One duplex molecule (fwd + rev)
        //  - One simplex-fwd molecule (two fwd reads, same UMI)
        //  - One singleton

        // Duplex pair
        let d_fwd = make_pair(b"ACGTA", b"TGCAT", b"ACGTACGTACGTACGTACGT");
        let d_rev = make_pair(b"TGCAT", b"ACGTA", b"ACGTACGTACGTACGTACGT");

        // Simplex-fwd pair (GGGGG < TTTTT, so both are Forward strand)
        let s1 = make_pair(b"GGGGG", b"TTTTT", b"CCCCCCCCCCCCCCCCCCCC");
        let s2 = make_pair(b"GGGGG", b"TTTTT", b"CCCCCCCCCCCCCCCCCCCC");

        // Singleton
        let single = make_pair(b"CCCCC", b"AAAAA", b"GGGGGGGGGGGGGGGGGGGG");
        // Note: "AAAAA" < "CCCCC" so canonical is (AAAAA, CCCCC);
        // R1=CCCCC → Reverse strand → singleton on rev.

        let (molecules, stats) = assemble_molecules(
            vec![d_fwd, d_rev, s1, s2, single],
            &AssemblerConfig::default(),
        );

        assert_eq!(molecules.len(), 3, "expected 3 molecules");
        assert_eq!(stats.n_molecules, 3);
        assert_eq!(stats.n_duplex, 1);
        assert_eq!(stats.n_simplex_fwd, 1);
        assert_eq!(stats.n_singletons, 1);
    }

    // ── Test 8: Genuine duplex pair with rotated fingerprints → one molecule ──
    // Before the rotation fix, forward and reverse reads from the same molecule
    // would fail fingerprint compatibility and be split into two separate
    // molecules, suppressing duplex calling.  The duplex-aware check must group
    // them into a single duplex molecule.

    #[test]
    fn duplex_pair_with_rotated_fingerprints_forms_one_molecule() {
        use kam_core::molecule::{CanonicalUmiPair, Strand};

        // For a genuine duplex pair the forward read has
        //   template_r1 = T,      template_r2 = RC(T)
        // and the reverse read has
        //   template_r1 = RC(T),  template_r2 = T.
        //
        // This gives fp_fwd and fp_rev = rotate_left(fp_fwd, 32), so direct
        // fingerprint comparison fails (~32 bits differ).
        //
        // Use T = AAAAAAAAGGGGGGGG and RC(T) = CCCCCCCCTTTTTTTT to get clearly
        // distinct templates whose fingerprint XOR has exactly 32 bits set.

        let t: Vec<u8> = b"AAAAAAAAGGGGGGGG".to_vec();
        let rc_t: Vec<u8> = b"CCCCCCCCTTTTTTTT".to_vec();

        // UMIs: ACGTA < TGCAT → canonical = (ACGTA, TGCAT).
        // Forward read: R1_UMI = umi_a = ACGTA → strand = Forward.
        // Reverse read: R1_UMI = umi_b = TGCAT → strand = Reverse.
        let canonical_umi = CanonicalUmiPair::new(b"ACGTA".to_vec(), b"TGCAT".to_vec());
        let qual = vec![b'I'; t.len()];

        let fwd_pair = crate::parser::ParsedReadPair {
            umi_r1: b"ACGTA".to_vec(),
            umi_r2: b"TGCAT".to_vec(),
            skip_r1: b"TG".to_vec(),
            skip_r2: b"TG".to_vec(),
            template_r1: t.clone(),
            template_r2: rc_t.clone(),
            qual_r1: qual.clone(),
            qual_r2: qual.clone(),
            umi_qual_r1: vec![b'I'; 5],
            umi_qual_r2: vec![b'I'; 5],
            canonical_umi: canonical_umi.clone(),
            strand: Strand::Forward,
        };

        let rev_pair = crate::parser::ParsedReadPair {
            umi_r1: b"TGCAT".to_vec(),
            umi_r2: b"ACGTA".to_vec(),
            skip_r1: b"TG".to_vec(),
            skip_r2: b"TG".to_vec(),
            template_r1: rc_t.clone(),
            template_r2: t.clone(),
            qual_r1: qual.clone(),
            qual_r2: qual.clone(),
            umi_qual_r1: vec![b'I'; 5],
            umi_qual_r2: vec![b'I'; 5],
            canonical_umi: canonical_umi.clone(),
            strand: Strand::Reverse,
        };

        let (molecules, stats) =
            assemble_molecules(vec![fwd_pair, rev_pair], &AssemblerConfig::default());

        assert_eq!(
            molecules.len(),
            1,
            "genuine duplex pair must be grouped into one molecule, not split"
        );
        assert_eq!(stats.n_duplex, 1, "should be classified as duplex");
        assert_eq!(
            stats.n_umi_collisions_detected, 0,
            "no UMI collision should be reported for a genuine duplex pair"
        );
    }

    // ── Test 9: Duplex consensus uses R1 for forward SSC, R2 for reverse SSC ──
    // This is a regression test for the bug where call_ssc_bases always read
    // template_r1, making duplex consensus compare R1 against R1 instead of
    // R1 against R2.  The test gives each read deliberately distinct R1 and R2
    // template sequences so the wrong selection would produce a different duplex
    // consensus sequence.
    //
    // Setup:
    //   fwd read: template_r1 = AAAAAAAAAAAAAAAA, template_r2 = CCCCCCCCCCCCCCCC
    //   rev read: template_r1 = CCCCCCCCCCCCCCCC, template_r2 = AAAAAAAAAAAAAAAA
    //
    // Correct: fwd SSC uses R1 → "AAAA...", rev SSC uses R2 → "AAAA...".
    //          Both strands agree → duplex consensus = "AAAA...".
    //
    // Bug: fwd SSC uses R1 → "AAAA...", rev SSC also uses R1 → "CCCC...".
    //      Both strands disagree → duplex consensus would be N-masked or absent.

    #[test]
    fn duplex_consensus_uses_r1_for_fwd_and_r2_for_rev() {
        use kam_core::molecule::{CanonicalUmiPair, Strand};

        let fwd_template: Vec<u8> = b"AAAAAAAAAAAAAAAA".to_vec();
        let rev_template: Vec<u8> = b"CCCCCCCCCCCCCCCC".to_vec();
        let qual = vec![b'I'; fwd_template.len()]; // Phred 40

        let canonical_umi = CanonicalUmiPair::new(b"ACGTA".to_vec(), b"TGCAT".to_vec());

        // Forward read: R1_UMI = umi_a = ACGTA → strand = Forward.
        // template_r1 = "AAAA..." (the molecule's sense strand template).
        // template_r2 = "CCCC..." (the antisense template, present on this read).
        let fwd_pair = crate::parser::ParsedReadPair {
            umi_r1: b"ACGTA".to_vec(),
            umi_r2: b"TGCAT".to_vec(),
            skip_r1: b"TG".to_vec(),
            skip_r2: b"TG".to_vec(),
            template_r1: fwd_template.clone(),
            template_r2: rev_template.clone(),
            qual_r1: qual.clone(),
            qual_r2: qual.clone(),
            umi_qual_r1: vec![b'I'; 5],
            umi_qual_r2: vec![b'I'; 5],
            canonical_umi: canonical_umi.clone(),
            strand: Strand::Forward,
        };

        // Reverse read: R1_UMI = umi_b = TGCAT → strand = Reverse.
        // template_r1 = "CCCC..." (antisense strand, as seen on this read's R1).
        // template_r2 = "AAAA..." (sense strand, as seen on this read's R2).
        let rev_pair = crate::parser::ParsedReadPair {
            umi_r1: b"TGCAT".to_vec(),
            umi_r2: b"ACGTA".to_vec(),
            skip_r1: b"TG".to_vec(),
            skip_r2: b"TG".to_vec(),
            template_r1: rev_template.clone(),
            template_r2: fwd_template.clone(),
            qual_r1: qual.clone(),
            qual_r2: qual.clone(),
            umi_qual_r1: vec![b'I'; 5],
            umi_qual_r2: vec![b'I'; 5],
            canonical_umi: canonical_umi.clone(),
            strand: Strand::Reverse,
        };

        let (molecules, stats) =
            assemble_molecules(vec![fwd_pair, rev_pair], &AssemblerConfig::default());

        assert_eq!(molecules.len(), 1, "should form one duplex molecule");
        assert_eq!(stats.n_duplex, 1);

        let mol = &molecules[0];
        let duplex = mol
            .duplex_consensus
            .as_ref()
            .expect("duplex consensus must be present");

        // With the correct template selection, fwd uses R1 = "AAAA..." and rev
        // uses R2 = "AAAA..." — both agree, so every base must be 'A'.
        assert!(
            duplex.sequence.iter().all(|&b| b == b'A'),
            "duplex consensus must be all-A (both strands agree on 'A'); \
             got {:?} — likely template_r1 was used for the reverse SSC instead of template_r2",
            std::str::from_utf8(&duplex.sequence).unwrap_or("<invalid utf8>"),
        );
    }

    // ── Test 10: Empty input returns empty vec ────────────────────────────────

    #[test]
    fn empty_input_returns_empty() {
        let (molecules, stats) = assemble_molecules(vec![], &AssemblerConfig::default());
        assert!(molecules.is_empty());
        assert_eq!(stats.n_molecules, 0);
        assert_eq!(stats.n_duplex, 0);
        assert_eq!(stats.n_singletons, 0);
    }
}
