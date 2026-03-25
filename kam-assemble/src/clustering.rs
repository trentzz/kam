//! Hamming-distance UMI pair clustering.
//!
//! Groups [`CanonicalUmiPair`]s whose combined UMI sequences differ by at most
//! `max_distance` positions. Uses a directional strategy: a higher-count UMI
//! absorbs a lower-count neighbour, but not the other way around. There is no
//! transitive chaining — if A absorbs B and B is near C, A will only absorb C
//! if A and C are themselves within `max_distance`.

use kam_core::molecule::CanonicalUmiPair;

/// Compute the Hamming distance between two canonical UMI pairs.
///
/// Compares `umi_a` of `a` against `umi_a` of `b`, and `umi_b` of `a` against
/// `umi_b` of `b`, then returns the total number of mismatching positions across
/// both arms. Positions beyond the shorter arm are not counted (zip stops at the
/// shorter of the two). UMIs of any length are supported.
///
/// # Example
/// ```
/// use kam_core::molecule::CanonicalUmiPair;
/// use kam_assemble::clustering::umi_pair_hamming_distance;
///
/// let a = CanonicalUmiPair::new(b"ACGTA".to_vec(), b"TGCAT".to_vec());
/// let b = CanonicalUmiPair::new(b"ACGTA".to_vec(), b"TGCAT".to_vec());
/// assert_eq!(umi_pair_hamming_distance(&a, &b), 0);
///
/// let c = CanonicalUmiPair::new(b"ACGTT".to_vec(), b"TGCAT".to_vec()); // last base of umi_a differs
/// assert_eq!(umi_pair_hamming_distance(&a, &c), 1);
/// ```
pub fn umi_pair_hamming_distance(a: &CanonicalUmiPair, b: &CanonicalUmiPair) -> u32 {
    let dist_a = a
        .umi_a
        .iter()
        .zip(b.umi_a.iter())
        .filter(|(x, y)| x != y)
        .count() as u32;
    let dist_b = a
        .umi_b
        .iter()
        .zip(b.umi_b.iter())
        .filter(|(x, y)| x != y)
        .count() as u32;
    dist_a + dist_b
}

/// Group canonical UMI pairs by Hamming distance using directional clustering.
///
/// **Algorithm**: pairs are sorted by read count (descending). Each pair is
/// either absorbed into an existing cluster whose seed is within `max_distance`,
/// or it starts a new cluster. Once a pair has been absorbed it cannot act as a
/// seed for another cluster. The comparison is always between the candidate and
/// the *seed* (first element) of an existing cluster — there is no transitive
/// chaining.
///
/// Returns a `Vec` of groups where each group is a `Vec<usize>` of indices into
/// the original `pairs` slice. Groups are returned in seed-first order, with
/// seeds having the highest count appearing earliest.
///
/// # Arguments
/// * `pairs` — slice of `(CanonicalUmiPair, read_count)` tuples.
/// * `max_distance` — maximum Hamming distance for two pairs to be merged.
///   Use `0` to skip clustering (each pair becomes its own singleton group).
///
/// # Example
/// ```
/// use kam_core::molecule::CanonicalUmiPair;
/// use kam_assemble::clustering::cluster_umi_pairs;
///
/// let pairs = vec![
///     (CanonicalUmiPair::new(b"ACGTA".to_vec(), b"TGCAT".to_vec()), 10_u32),
///     (CanonicalUmiPair::new(b"ACGTT".to_vec(), b"TGCAT".to_vec()), 2_u32),  // 1 mismatch from idx 0
/// ];
/// let groups = cluster_umi_pairs(&pairs, 1);
/// assert_eq!(groups.len(), 1);
/// assert_eq!(groups[0].len(), 2);
/// ```
pub fn cluster_umi_pairs(pairs: &[(CanonicalUmiPair, u32)], max_distance: u32) -> Vec<Vec<usize>> {
    // Build an index sorted by count descending (stable sort preserves input
    // order for ties, making output deterministic).
    let mut order: Vec<usize> = (0..pairs.len()).collect();
    order.sort_by(|&i, &j| pairs[j].1.cmp(&pairs[i].1));

    // `cluster_of[i]` holds the cluster index assigned to input index i, or
    // None if it has not yet been assigned.
    let mut cluster_of: Vec<Option<usize>> = vec![None; pairs.len()];
    // Each entry is (seed_input_index, members).
    let mut clusters: Vec<(usize, Vec<usize>)> = Vec::new();

    for &idx in &order {
        if cluster_of[idx].is_some() {
            // Already absorbed — skip.
            continue;
        }

        // Try to find an existing cluster whose seed is within max_distance.
        let mut found = None;
        for (cluster_idx, (seed_idx, _members)) in clusters.iter().enumerate() {
            let dist = umi_pair_hamming_distance(&pairs[*seed_idx].0, &pairs[idx].0);
            if dist <= max_distance {
                found = Some(cluster_idx);
                break;
            }
        }

        match found {
            Some(cluster_idx) => {
                cluster_of[idx] = Some(cluster_idx);
                clusters[cluster_idx].1.push(idx);
            }
            None => {
                // Start a new cluster with this pair as the seed.
                let new_cluster_idx = clusters.len();
                cluster_of[idx] = Some(new_cluster_idx);
                clusters.push((idx, vec![idx]));
            }
        }
    }

    clusters
        .into_iter()
        .map(|(_seed, members)| members)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pair(a: &[u8], b: &[u8]) -> CanonicalUmiPair {
        CanonicalUmiPair::new(a.to_vec(), b.to_vec())
    }

    // Test 1: Hamming distance of identical pairs is 0.
    #[test]
    fn hamming_identical_is_zero() {
        let a = pair(b"ACGTA", b"TGCAT");
        assert_eq!(umi_pair_hamming_distance(&a, &a), 0);
    }

    // Test 2: Hamming distance of pairs differing in 1 base of umi_a is 1.
    #[test]
    fn hamming_one_mismatch_in_umi_a() {
        let a = pair(b"ACGTA", b"TGCAT");
        // Last base of umi_a changes: A -> T.  umi_b stays the same.
        // After canonicalisation: umi_a="ACGTT", umi_b="TGCAT"
        let b = pair(b"ACGTT", b"TGCAT");
        assert_eq!(umi_pair_hamming_distance(&a, &b), 1);
    }

    // Test 3: Hamming distance of pairs differing in 1 base of each UMI is 2.
    #[test]
    fn hamming_one_mismatch_in_each_umi() {
        let a = pair(b"ACGTA", b"TGCAT");
        // 1 mismatch in umi_a (last base) + 1 mismatch in umi_b (last base).
        // "ACGTT" < "TGCAG" so canonical order unchanged.
        let b = pair(b"ACGTT", b"TGCAG");
        assert_eq!(umi_pair_hamming_distance(&a, &b), 2);
    }

    // Test 4: Clustering with max_distance=0 produces one group per unique pair.
    #[test]
    fn cluster_max_distance_zero_one_group_per_pair() {
        let pairs = vec![
            (pair(b"ACGTA", b"TGCAT"), 10_u32),
            (pair(b"ACGTT", b"TGCAT"), 5_u32),
            (pair(b"AAAAA", b"TTTTT"), 3_u32),
        ];
        let groups = cluster_umi_pairs(&pairs, 0);
        assert_eq!(groups.len(), 3, "each distinct pair is its own group");
        // Every group has exactly one member.
        assert!(groups.iter().all(|g| g.len() == 1));
    }

    // Test 5: Clustering merges pairs within distance 1.
    #[test]
    fn cluster_merges_within_distance_one() {
        let pairs = vec![
            (pair(b"ACGTA", b"TGCAT"), 10_u32),
            (pair(b"ACGTT", b"TGCAT"), 2_u32), // 1 mismatch from index 0
        ];
        let groups = cluster_umi_pairs(&pairs, 1);
        assert_eq!(groups.len(), 1, "the two pairs should form one cluster");
        assert_eq!(groups[0].len(), 2);
        // Both original indices are present.
        let mut members = groups[0].clone();
        members.sort_unstable();
        assert_eq!(members, vec![0, 1]);
    }

    // Test 6: Directional — high-count seed absorbs low-count neighbour, not vice versa.
    #[test]
    fn cluster_directional_high_absorbs_low() {
        // Index 0: high-count; index 1: low-count, 1 mismatch away.
        // Index 2: distinct pair that is 1 mismatch from index 1 only.
        //
        // With directional clustering the processing order is 0, 1, 2 (by
        // count desc).  Index 1 gets absorbed into cluster seeded by 0.
        // Index 2's seed is index 0; distance(0, 2) must exceed max_distance
        // for it to remain in its own cluster.
        //
        // We make index 2 only close to index 1, not to index 0.
        let high = pair(b"ACGTA", b"TGCAT"); // seed
        let low = pair(b"ACGTT", b"TGCAT"); // 1 mismatch from high  → absorbed by high
        let other = pair(b"AAAAA", b"CCCCC"); // far from high, far from low

        let pairs = vec![(high.clone(), 20_u32), (low, 3_u32), (other, 1_u32)];
        let groups = cluster_umi_pairs(&pairs, 1);

        // Cluster with seed 0 (high) should contain indices 0 and 1.
        // Cluster with seed 2 (other) should contain only index 2.
        assert_eq!(groups.len(), 2);

        // Find the group that contains index 0.
        let group_of_high = groups
            .iter()
            .find(|g| g.contains(&0))
            .expect("high must be in a group");
        assert!(group_of_high.contains(&1), "low should be absorbed by high");
        assert!(
            !group_of_high.contains(&2),
            "other should not be in high's group"
        );

        let group_of_other = groups
            .iter()
            .find(|g| g.contains(&2))
            .expect("other must be in a group");
        assert_eq!(group_of_other.len(), 1);
    }

    // Test 7: No transitive chaining — A↔B within distance, B↔C within distance,
    // but A↔C NOT within distance: A and C must be in separate clusters.
    #[test]
    fn cluster_no_transitive_chaining() {
        // A: "AAAAA+TTTTT"  (seed, highest count)
        // B: "AAAAC+TTTTT"  (1 mismatch from A in umi_a)
        // C: "AAACC+TTTTT"  (1 mismatch from B in umi_a, but 2 mismatches from A)
        //
        // With max_distance=1: A absorbs B; C is checked against A's seed only,
        // distance(A,C)=2 > 1, so C starts its own cluster.
        let a = pair(b"AAAAA", b"TTTTT");
        let b_pair = pair(b"AAAAC", b"TTTTT"); // distance(A,B)=1
        let c = pair(b"AAACC", b"TTTTT"); // distance(A,C)=2, distance(B,C)=1

        let pairs = vec![(a, 100_u32), (b_pair, 10_u32), (c, 5_u32)];
        let groups = cluster_umi_pairs(&pairs, 1);

        assert_eq!(groups.len(), 2, "A+B in one cluster, C alone");

        let group_of_a = groups
            .iter()
            .find(|g| g.contains(&0))
            .expect("A must be in a group");
        assert!(
            group_of_a.contains(&1),
            "B should be absorbed into A's cluster"
        );
        assert!(
            !group_of_a.contains(&2),
            "C must NOT be absorbed transitively"
        );

        let group_of_c = groups
            .iter()
            .find(|g| g.contains(&2))
            .expect("C must be in a group");
        assert_eq!(group_of_c.len(), 1);
    }

    // Test 8: Single-element input returns a single group.
    #[test]
    fn cluster_single_element_returns_single_group() {
        let pairs = vec![(pair(b"ACGTA", b"TGCAT"), 5_u32)];
        let groups = cluster_umi_pairs(&pairs, 1);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0], vec![0]);
    }

    // Test 9: Hamming distance with 12 bp UMIs — identical pairs give distance 0.
    #[test]
    fn hamming_twelve_bp_umi_identical() {
        let a = pair(b"ACGTACGTACGT", b"TGCATGCATGCA");
        assert_eq!(umi_pair_hamming_distance(&a, &a), 0);
    }

    // Test 10: Hamming distance with 12 bp UMIs — two mismatches detected.
    #[test]
    fn hamming_twelve_bp_umi_two_mismatches() {
        // "ACGTACGTACGT" vs "ACGTACGTACCC" — last 2 bases of umi_a differ (T→C, T→C = 2 mismatches).
        // "TGCATGCA" < "TGCATGCATGCA" — using distinct long umis to ensure canonical ordering is stable.
        let a = pair(b"ACGTACGTACGT", b"TGCATGCATGCA");
        // Two mismatches in umi_a.
        let b = pair(b"ACGTACGTACCC", b"TGCATGCATGCA");
        assert_eq!(umi_pair_hamming_distance(&a, &b), 2);
    }

    // Test 11: Clustering with 12 bp UMIs — pairs within Hamming-1 merge.
    #[test]
    fn cluster_twelve_bp_umi_merges_within_distance_one() {
        // Two pairs: last base of umi_a differs by 1 (T → A = 1 mismatch).
        let pairs = vec![
            (pair(b"ACGTACGTACGT", b"TGCATGCATGCA"), 10_u32),
            (pair(b"ACGTACGTACGA", b"TGCATGCATGCA"), 2_u32), // 1 mismatch in umi_a
        ];
        let groups = cluster_umi_pairs(&pairs, 1);
        assert_eq!(groups.len(), 1, "Hamming-1 12 bp UMI pairs should merge");
        assert_eq!(groups[0].len(), 2);
    }

    // Test 12: Hamming distance with 9 bp UMIs — one mismatch detected.
    #[test]
    fn hamming_nine_bp_umi_one_mismatch() {
        let a = pair(b"AAAAAAAAA", b"TTTTTTTTT");
        // Last base of umi_a changes: A → C = 1 mismatch.
        let b = pair(b"AAAAAAAAC", b"TTTTTTTTT");
        assert_eq!(umi_pair_hamming_distance(&a, &b), 1);
    }

    // Test 13: Clustering with 9 bp UMIs — distant pairs stay separate.
    #[test]
    fn cluster_nine_bp_umi_distinct_pairs_stay_separate() {
        let pairs = vec![
            (pair(b"AAAAAAAAA", b"TTTTTTTTT"), 10_u32),
            (pair(b"CCCCCCCCC", b"GGGGGGGGG"), 10_u32), // very distant
        ];
        let groups = cluster_umi_pairs(&pairs, 1);
        assert_eq!(groups.len(), 2, "distant 9 bp UMI pairs must not be merged");
    }
}
