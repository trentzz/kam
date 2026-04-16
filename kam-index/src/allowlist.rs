//! Allowlist-based k-mer filtering for panel-targeted indexing.
//!
//! Building an allowlist from target sequences bounds memory usage to the panel
//! size rather than sequencing depth: only k-mers that could plausibly arise
//! from the targeted regions are stored in the index.

use std::collections::HashSet;

use kam_core::kmer::{KmerIndex, MoleculeEvidence};

use crate::encode::{canonical, KmerIterator};

/// Build an allowlist of canonical k-mers from target sequences.
///
/// Every k-mer in every target sequence is canonicalised and added to the
/// returned set. The caller is responsible for providing sequences with any
/// desired flanking context already included — this function does not pad.
///
/// # Arguments
///
/// * `targets` — one or more target sequences (ASCII DNA bytes).
/// * `k` — k-mer length (must be in `1..=31`).
///
/// # Examples
///
/// ```
/// use kam_index::allowlist::build_allowlist;
///
/// let targets: &[&[u8]] = &[b"ACGTACGT"];
/// let al = build_allowlist(targets, 4);
/// // 5 k-mers in an 8-base sequence, but some may collapse after canonicalisation.
/// assert!(!al.is_empty());
/// ```
pub fn build_allowlist(targets: &[&[u8]], k: usize) -> HashSet<u64> {
    let mut set = HashSet::new();
    for &target in targets {
        for (_pos, kmer) in KmerIterator::new(target, k) {
            set.insert(canonical(kmer, k));
        }
    }
    set
}

/// A k-mer index that only stores k-mers present in an allowlist.
///
/// Wraps any inner [`KmerIndex`] implementation and silently discards
/// insertions for k-mers not in the allowlist. All reads delegate to the
/// inner index, so only allowed k-mers are ever returned.
///
/// This is the primary memory-bounding mechanism for targeted panel work:
/// memory scales with the number of panel k-mers, not sequencing depth.
///
/// # Example
///
/// ```
/// use std::collections::HashSet;
/// use kam_index::allowlist::{build_allowlist, FilteredKmerIndex};
/// use kam_index::HashKmerIndex;
/// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
///
/// let targets: &[&[u8]] = &[b"ACGTACGT"];
/// let al = build_allowlist(targets, 4);
/// let mut index = FilteredKmerIndex::new(HashKmerIndex::new(), al);
///
/// // Insert an evidence record for an allowed k-mer.
/// let first_kmer = build_allowlist(&[b"ACGT"], 4).into_iter().next().unwrap();
/// index.insert(first_kmer, MoleculeEvidence { n_molecules: 1, ..Default::default() });
/// assert_eq!(index.observed_count(), 1);
/// ```
#[derive(Debug, Clone)]
pub struct FilteredKmerIndex<T: KmerIndex> {
    inner: T,
    allowlist: HashSet<u64>,
}

impl<T: KmerIndex> FilteredKmerIndex<T> {
    /// Create a new [`FilteredKmerIndex`] wrapping `inner` with the given `allowlist`.
    ///
    /// # Example
    ///
    /// ```
    /// use std::collections::HashSet;
    /// use kam_index::allowlist::FilteredKmerIndex;
    /// use kam_index::HashKmerIndex;
    ///
    /// let al: HashSet<u64> = [1u64, 2, 3].into_iter().collect();
    /// let index = FilteredKmerIndex::new(HashKmerIndex::new(), al);
    /// assert_eq!(index.allowlist_size(), 3);
    /// ```
    pub fn new(inner: T, allowlist: HashSet<u64>) -> Self {
        Self { inner, allowlist }
    }

    /// Return the number of k-mers in the allowlist.
    ///
    /// # Example
    ///
    /// ```
    /// use std::collections::HashSet;
    /// use kam_index::allowlist::FilteredKmerIndex;
    /// use kam_index::HashKmerIndex;
    ///
    /// let al: HashSet<u64> = [10u64, 20].into_iter().collect();
    /// let index = FilteredKmerIndex::new(HashKmerIndex::new(), al);
    /// assert_eq!(index.allowlist_size(), 2);
    /// ```
    pub fn allowlist_size(&self) -> usize {
        self.allowlist.len()
    }

    /// Return the number of allowed k-mers that have been observed (have evidence).
    ///
    /// # Example
    ///
    /// ```
    /// use std::collections::HashSet;
    /// use kam_index::allowlist::FilteredKmerIndex;
    /// use kam_index::HashKmerIndex;
    /// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
    ///
    /// let al: HashSet<u64> = [42u64].into_iter().collect();
    /// let mut index = FilteredKmerIndex::new(HashKmerIndex::new(), al);
    /// assert_eq!(index.observed_count(), 0);
    /// index.insert(42, MoleculeEvidence::default());
    /// assert_eq!(index.observed_count(), 1);
    /// ```
    pub fn observed_count(&self) -> usize {
        self.allowlist
            .iter()
            .filter(|&&k| self.inner.contains(k))
            .count()
    }

    /// Return `true` if `kmer` is in the allowlist, regardless of whether it
    /// has been observed.
    ///
    /// # Example
    ///
    /// ```
    /// use std::collections::HashSet;
    /// use kam_index::allowlist::FilteredKmerIndex;
    /// use kam_index::HashKmerIndex;
    ///
    /// let al: HashSet<u64> = [7u64].into_iter().collect();
    /// let index = FilteredKmerIndex::new(HashKmerIndex::new(), al);
    /// assert!(index.is_allowed(7));
    /// assert!(!index.is_allowed(8));
    /// ```
    pub fn is_allowed(&self, kmer: u64) -> bool {
        self.allowlist.contains(&kmer)
    }
}

impl<T: KmerIndex> KmerIndex for FilteredKmerIndex<T> {
    /// Return a reference to the evidence for `kmer`, or `None` if absent.
    fn get(&self, kmer: u64) -> Option<&MoleculeEvidence> {
        self.inner.get(kmer)
    }

    /// Insert evidence for `kmer`.
    ///
    /// This is a **no-op** if `kmer` is not in the allowlist — the insertion
    /// is silently dropped without an error or panic.
    fn insert(&mut self, kmer: u64, evidence: MoleculeEvidence) {
        if self.allowlist.contains(&kmer) {
            self.inner.insert(kmer, evidence);
        }
    }

    /// Return `true` if the index contains any evidence for `kmer`.
    fn contains(&self, kmer: u64) -> bool {
        self.inner.contains(kmer)
    }

    /// Return `n_molecules` for `kmer`, or `0` if absent.
    fn molecule_count(&self, kmer: u64) -> u32 {
        self.inner.molecule_count(kmer)
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode::encode_kmer;
    use crate::hash_index::HashKmerIndex;
    use kam_core::kmer::KmerIndex;

    fn ev(n: u32) -> MoleculeEvidence {
        MoleculeEvidence {
            n_molecules: n,
            ..Default::default()
        }
    }

    // Test 1: build_allowlist from a simple target produces expected k-mers.
    #[test]
    fn build_allowlist_simple_target() {
        let target: &[u8] = b"ACGTACGT";
        let k = 4;
        let al = build_allowlist(&[target], k);

        // An 8-base sequence has 5 k-mers; after canonicalisation some may collapse.
        assert!(!al.is_empty());
        assert!(al.len() <= 5);

        // Every k-mer extracted from the sequence should be in the allowlist.
        for (_pos, kmer) in KmerIterator::new(target, k) {
            assert!(al.contains(&canonical(kmer, k)));
        }
    }

    // Test 2: build_allowlist canonicalises — forward and revcomp map to same entry.
    #[test]
    fn build_allowlist_canonicalises_kmers() {
        // "AAAC" and "GTTT" are reverse complements of each other.
        let k = 4;
        let fwd = encode_kmer(b"AAAC").unwrap();
        let rev = encode_kmer(b"GTTT").unwrap();

        let al_fwd = build_allowlist(&[b"AAAC"], k);
        let al_rev = build_allowlist(&[b"GTTT"], k);

        // Both should produce the same canonical form.
        let can_fwd = canonical(fwd, k);
        let can_rev = canonical(rev, k);
        assert_eq!(can_fwd, can_rev);
        assert!(al_fwd.contains(&can_fwd));
        assert!(al_rev.contains(&can_rev));
        // The two allowlists contain the same single element.
        assert_eq!(al_fwd, al_rev);
    }

    // Test 3: FilteredKmerIndex insert + get works for allowed k-mers.
    #[test]
    fn filtered_insert_and_get_allowed_kmer() {
        let al: HashSet<u64> = [100u64, 200].into_iter().collect();
        let mut idx = FilteredKmerIndex::new(HashKmerIndex::new(), al);

        idx.insert(100, ev(3));
        let got = idx
            .get(100)
            .expect("should have evidence for allowed k-mer");
        assert_eq!(got.n_molecules, 3);
    }

    // Test 4: FilteredKmerIndex insert is no-op for non-allowed k-mers.
    #[test]
    fn filtered_insert_noop_for_disallowed_kmer() {
        let al: HashSet<u64> = [100u64].into_iter().collect();
        let mut idx = FilteredKmerIndex::new(HashKmerIndex::new(), al);

        idx.insert(999, ev(5)); // 999 is not in the allowlist
        assert!(!idx.contains(999));
        assert!(idx.get(999).is_none());
    }

    // Test 5: allowlist_size returns correct count.
    #[test]
    fn allowlist_size_correct() {
        let al: HashSet<u64> = [1u64, 2, 3, 4, 5].into_iter().collect();
        let idx = FilteredKmerIndex::new(HashKmerIndex::new(), al);
        assert_eq!(idx.allowlist_size(), 5);
    }

    // Test 6: is_allowed true for allowed, false for others.
    #[test]
    fn is_allowed_true_and_false() {
        let al: HashSet<u64> = [42u64, 43].into_iter().collect();
        let idx = FilteredKmerIndex::new(HashKmerIndex::new(), al);
        assert!(idx.is_allowed(42));
        assert!(idx.is_allowed(43));
        assert!(!idx.is_allowed(44));
        assert!(!idx.is_allowed(0));
    }

    // Test 7: Multiple targets: allowlist is union of all target k-mers.
    #[test]
    fn build_allowlist_union_of_multiple_targets() {
        let k = 4;
        let al_a = build_allowlist(&[b"AAAA"], k);
        let al_b = build_allowlist(&[b"CCCC"], k);
        let al_both = build_allowlist(&[b"AAAA", b"CCCC"], k);

        // Union must contain everything from both individual allowlists.
        for &kmer in &al_a {
            assert!(al_both.contains(&kmer));
        }
        for &kmer in &al_b {
            assert!(al_both.contains(&kmer));
        }
        // Both single-base sequences each produce 1 canonical k-mer.
        assert_eq!(al_both.len(), al_a.len() + al_b.len());
    }

    // Test 8: observed_count starts at 0, increases with inserts.
    #[test]
    fn observed_count_starts_zero_and_increments() {
        let al: HashSet<u64> = [10u64, 20, 30].into_iter().collect();
        let mut idx = FilteredKmerIndex::new(HashKmerIndex::new(), al);

        assert_eq!(idx.observed_count(), 0);

        idx.insert(10, ev(1));
        assert_eq!(idx.observed_count(), 1);

        idx.insert(20, ev(2));
        assert_eq!(idx.observed_count(), 2);

        // Re-inserting the same k-mer doesn't increase observed_count.
        idx.insert(10, ev(99));
        assert_eq!(idx.observed_count(), 2);

        // Inserting a non-allowed k-mer doesn't change observed_count.
        idx.insert(999, ev(1));
        assert_eq!(idx.observed_count(), 2);
    }

    // ── New edge-case tests ──────────────────────────────────────────────────

    // Test 9: Sequence shorter than k produces no k-mers (no panic).
    #[test]
    fn build_allowlist_seq_shorter_than_k_produces_empty() {
        let al = build_allowlist(&[b"AC"], 4);
        assert!(al.is_empty(), "a 2-base sequence cannot yield any 4-mers");
    }

    // Test 10: Sequence exactly length k produces exactly one canonical k-mer.
    #[test]
    fn build_allowlist_seq_exactly_k_produces_one_kmer() {
        let al = build_allowlist(&[b"ACGT"], 4);
        assert_eq!(al.len(), 1, "a 4-base sequence yields exactly one 4-mer");
        // ACGT is a palindrome (its own reverse complement), so the canonical
        // form equals the raw encoding.
        let expected = canonical(encode_kmer(b"ACGT").expect("valid k-mer"), 4);
        assert!(al.contains(&expected));
    }

    // Test 11: Sequence of all N produces no k-mers.
    // KmerIterator skips N bases (encode_base returns None), so no valid
    // k-mer window ever forms.
    #[test]
    fn build_allowlist_all_n_produces_empty() {
        let al = build_allowlist(&[b"NNNNNNNN"], 4);
        assert!(
            al.is_empty(),
            "all-N sequence should produce no valid k-mers"
        );
    }

    // Test 12: Palindromic k-mer (its own reverse complement) appears once.
    // ACGT is self-complementary. Canonicalising it yields the same value,
    // so the set contains exactly one entry, not two.
    #[test]
    fn build_allowlist_palindromic_kmer_appears_once() {
        let al = build_allowlist(&[b"ACGT"], 4);
        assert_eq!(al.len(), 1, "palindromic k-mer should appear once in set");

        let enc = encode_kmer(b"ACGT").expect("valid k-mer");
        let rc = crate::encode::reverse_complement(enc, 4);
        assert_eq!(enc, rc, "ACGT must be its own reverse complement");
    }

    // Test 13: Same sequence provided twice gives the same allowlist size.
    // build_allowlist collects into a HashSet, so duplicates collapse.
    #[test]
    fn build_allowlist_duplicate_sequence_same_size() {
        let single = build_allowlist(&[b"ACGTACGT"], 4);
        let doubled = build_allowlist(&[b"ACGTACGT", b"ACGTACGT"], 4);
        assert_eq!(
            single.len(),
            doubled.len(),
            "duplicate target sequences must not inflate the allowlist"
        );
    }

    // Test 14: Two distinct sequences with overlapping k-mers. The union
    // contains only the distinct canonical k-mers across both.
    #[test]
    fn build_allowlist_overlapping_kmers_correct_union_size() {
        let k = 3;
        // "ACGT" yields k-mers: ACG, CGT (2 raw k-mers)
        // "CGTC" yields k-mers: CGT, GTC (2 raw k-mers)
        // CGT is shared. After canonicalisation the union should have 3
        // distinct entries (ACG/CGT overlap).
        let al_a = build_allowlist(&[b"ACGT"], k);
        let al_b = build_allowlist(&[b"CGTC"], k);
        let al_both = build_allowlist(&[b"ACGT", b"CGTC"], k);

        // The union must be at most the sum of the individual sizes but also
        // at least the larger of the two.
        assert!(al_both.len() <= al_a.len() + al_b.len());
        assert!(al_both.len() >= al_a.len().max(al_b.len()));

        // Every k-mer from each individual list is in the union.
        for &km in &al_a {
            assert!(al_both.contains(&km));
        }
        for &km in &al_b {
            assert!(al_both.contains(&km));
        }
    }

    // Test 15: k=1 edge case. Single-base k-mers. There are only 4 possible
    // raw 1-mers (A, C, G, T), and canonicalisation collapses A/T and C/G.
    #[test]
    fn build_allowlist_k_equals_one() {
        // "ACGT" has 4 raw 1-mers. After canonicalisation:
        //   A (0b00) and T (0b11) → canonical = min(A, T) = A
        //   C (0b01) and G (0b10) → canonical = min(C, G) = C
        // So only 2 distinct canonical 1-mers.
        let al = build_allowlist(&[b"ACGT"], 1);
        assert_eq!(al.len(), 2, "ACGT with k=1 should yield 2 canonical 1-mers");
    }

    // Test 16: Very long sequence (1000 bp) produces (len - k + 1) raw k-mers.
    // After canonicalisation, the set size should be <= (len - k + 1).
    #[test]
    fn build_allowlist_long_sequence_correct_upper_bound() {
        let k = 5;
        // Build a 1000-base sequence from a repeating pattern that avoids N.
        let pattern = b"ACGTACGTACGTACGTACGT"; // 20 bases
        let long_seq: Vec<u8> = pattern.iter().cycle().take(1000).copied().collect();
        let al = build_allowlist(&[&long_seq], k);

        let max_raw_kmers = 1000 - k + 1;
        assert!(
            !al.is_empty(),
            "1000-base sequence must produce some k-mers"
        );
        assert!(
            al.len() <= max_raw_kmers,
            "allowlist cannot exceed number of raw k-mers ({max_raw_kmers})"
        );
    }

    // Test 17: FilteredKmerIndex molecule_count returns 0 for disallowed k-mer
    // even after attempted insert.
    #[test]
    fn filtered_molecule_count_zero_for_disallowed() {
        let al: HashSet<u64> = [10u64].into_iter().collect();
        let mut idx = FilteredKmerIndex::new(HashKmerIndex::new(), al);
        idx.insert(999, ev(42));
        assert_eq!(
            idx.molecule_count(999),
            0,
            "disallowed k-mer must report zero molecules"
        );
    }

    // Test 18: Empty targets slice produces an empty allowlist (no panic).
    #[test]
    fn build_allowlist_empty_targets_produces_empty() {
        let targets: &[&[u8]] = &[];
        let al = build_allowlist(targets, 4);
        assert!(
            al.is_empty(),
            "no target sequences should yield an empty allowlist"
        );
    }

    // Test 19: A single N in the middle of a sequence splits the k-mer window.
    // "ACG N ACG" with k=3 yields ACG from both halves but only one canonical
    // entry because they are the same k-mer.
    #[test]
    fn build_allowlist_n_in_middle_splits_window() {
        let k = 3;
        // "ACGNACG" has 7 bases. Without the N we would get 5 k-mers.
        // With the N splitting at position 3, we get:
        //   Left half "ACG" → 1 k-mer (ACG)
        //   Right half "ACG" → 1 k-mer (ACG)
        // Both are the same, so the allowlist has at most 1 canonical entry.
        let al = build_allowlist(&[b"ACGNACG"], k);
        assert!(
            !al.is_empty(),
            "valid bases on both sides of N should produce k-mers"
        );
        // The allowlist should have only the canonical form of ACG (and its RC CGT
        // collapses into whichever is smaller).
        assert_eq!(
            al.len(),
            1,
            "both halves yield the same k-mer, so the set has one entry"
        );
    }

    // Test 20: Lowercase bases are handled correctly by build_allowlist.
    // encode_base accepts both 'a' and 'A', so lowercase input should work.
    #[test]
    fn build_allowlist_lowercase_bases_accepted() {
        let k = 4;
        let upper = build_allowlist(&[b"ACGT"], k);
        let lower = build_allowlist(&[b"acgt"], k);
        assert_eq!(
            upper, lower,
            "uppercase and lowercase input must produce identical allowlists"
        );
    }

    // Test 21: FilteredKmerIndex with an empty allowlist rejects all inserts.
    #[test]
    fn filtered_empty_allowlist_rejects_everything() {
        let al: HashSet<u64> = HashSet::new();
        let mut idx = FilteredKmerIndex::new(HashKmerIndex::new(), al);
        assert_eq!(idx.allowlist_size(), 0);

        idx.insert(1, ev(10));
        idx.insert(2, ev(20));
        assert_eq!(idx.observed_count(), 0);
        assert!(!idx.contains(1));
        assert!(!idx.contains(2));
    }

    // Test 22: build_allowlist with k=31 (maximum supported k-mer length).
    // A 31-base sequence should yield exactly one k-mer.
    #[test]
    fn build_allowlist_k_equals_31_max_supported() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACG"; // 31 bases
        assert_eq!(seq.len(), 31);
        let al = build_allowlist(&[seq], 31);
        assert_eq!(
            al.len(),
            1,
            "a 31-base sequence with k=31 yields exactly one canonical k-mer"
        );
    }

    // Test 23: Sequence with trailing N. Only the valid prefix contributes.
    #[test]
    fn build_allowlist_trailing_n_only_prefix_contributes() {
        let k = 3;
        // "ACGTN" has the N at position 4. Valid k-mers come from positions 0-2:
        // ACG (pos 0), CGT (pos 1). That is 2 raw k-mers.
        let al = build_allowlist(&[b"ACGTN"], k);
        assert!(!al.is_empty());
        // At most 2 canonical k-mers.
        assert!(al.len() <= 2);
    }
}
