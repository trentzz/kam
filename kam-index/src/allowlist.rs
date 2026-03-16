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
        self.allowlist.iter().filter(|&&k| self.inner.contains(k)).count()
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
        MoleculeEvidence { n_molecules: n, ..Default::default() }
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
        let got = idx.get(100).expect("should have evidence for allowed k-mer");
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
}
