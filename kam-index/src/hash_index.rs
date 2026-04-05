//! In-memory k-mer index backed by a [`HashMap`].
//!
//! [`HashKmerIndex`] is the primary index for targeted panel work where the
//! full k-mer space fits comfortably in RAM (~48 bytes per entry).

use std::collections::HashMap;

use kam_core::kmer::{KmerIndex, MoleculeEvidence};

/// In-memory k-mer index backed by [`HashMap`].
///
/// Stores one [`MoleculeEvidence`] record per k-mer (encoded as `u64`).
/// Fast insertion and lookup with O(1) average complexity.
/// Approximately 48 bytes per entry — best for targeted panel work.
///
/// # Example
/// ```
/// use kam_index::HashKmerIndex;
/// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
///
/// let mut index = HashKmerIndex::new();
/// let ev = MoleculeEvidence { n_molecules: 2, ..Default::default() };
/// index.insert(42, ev);
/// assert_eq!(index.molecule_count(42), 2);
/// assert!(index.contains(42));
/// ```
#[derive(Debug, Clone, Default)]
pub struct HashKmerIndex {
    map: HashMap<u64, MoleculeEvidence>,
}

impl HashKmerIndex {
    /// Create a new, empty [`HashKmerIndex`].
    ///
    /// # Example
    /// ```
    /// use kam_index::HashKmerIndex;
    ///
    /// let index = HashKmerIndex::new();
    /// assert!(index.is_empty());
    /// ```
    pub fn new() -> Self {
        Self {
            map: HashMap::new(),
        }
    }

    /// Create a new [`HashKmerIndex`] with pre-allocated capacity.
    ///
    /// Pre-allocation avoids rehashing when the expected number of k-mers is
    /// known up front. The behaviour is identical to [`HashKmerIndex::new`];
    /// only memory allocation is affected.
    ///
    /// # Example
    /// ```
    /// use kam_index::HashKmerIndex;
    ///
    /// let index = HashKmerIndex::with_capacity(1_000_000);
    /// assert!(index.is_empty());
    /// ```
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            map: HashMap::with_capacity(capacity),
        }
    }

    /// Return the number of k-mers stored in the index.
    ///
    /// # Example
    /// ```
    /// use kam_index::HashKmerIndex;
    /// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
    ///
    /// let mut index = HashKmerIndex::new();
    /// assert_eq!(index.len(), 0);
    /// index.insert(1, MoleculeEvidence::default());
    /// assert_eq!(index.len(), 1);
    /// ```
    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// Return `true` if the index contains no k-mers.
    ///
    /// # Example
    /// ```
    /// use kam_index::HashKmerIndex;
    ///
    /// let index = HashKmerIndex::new();
    /// assert!(index.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// Return a mutable reference to the evidence for `kmer`, inserting a
    /// default [`MoleculeEvidence`] if the k-mer is not yet present.
    ///
    /// This mirrors [`HashMap::entry`] and is useful for accumulating evidence
    /// across multiple observations without a separate lookup.
    ///
    /// # Example
    /// ```
    /// use kam_index::HashKmerIndex;
    /// use kam_core::kmer::MoleculeEvidence;
    ///
    /// let mut index = HashKmerIndex::new();
    /// let ev = index.entry(7);
    /// ev.n_molecules += 1;
    /// assert_eq!(index.entry(7).n_molecules, 1);
    /// ```
    pub fn entry(&mut self, kmer: u64) -> &mut MoleculeEvidence {
        self.map.entry(kmer).or_default()
    }

    /// Iterate over all k-mers and their associated [`MoleculeEvidence`].
    ///
    /// Iteration order is unspecified (HashMap order).
    ///
    /// # Example
    /// ```
    /// use kam_index::HashKmerIndex;
    /// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
    ///
    /// let mut index = HashKmerIndex::new();
    /// index.insert(1, MoleculeEvidence::default());
    /// index.insert(2, MoleculeEvidence::default());
    /// assert_eq!(index.iter().count(), 2);
    /// ```
    pub fn iter(&self) -> impl Iterator<Item = (&u64, &MoleculeEvidence)> {
        self.map.iter()
    }
}

impl KmerIndex for HashKmerIndex {
    /// Return a reference to the evidence for `kmer`, or `None` if absent.
    ///
    /// # Example
    /// ```
    /// use kam_index::HashKmerIndex;
    /// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
    ///
    /// let mut index = HashKmerIndex::new();
    /// assert!(index.get(99).is_none());
    /// index.insert(99, MoleculeEvidence { n_molecules: 1, ..Default::default() });
    /// assert!(index.get(99).is_some());
    /// ```
    fn get(&self, kmer: u64) -> Option<&MoleculeEvidence> {
        self.map.get(&kmer)
    }

    /// Insert or overwrite evidence for `kmer`.
    ///
    /// # Example
    /// ```
    /// use kam_index::HashKmerIndex;
    /// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
    ///
    /// let mut index = HashKmerIndex::new();
    /// index.insert(5, MoleculeEvidence { n_molecules: 3, ..Default::default() });
    /// assert_eq!(index.molecule_count(5), 3);
    /// ```
    fn insert(&mut self, kmer: u64, evidence: MoleculeEvidence) {
        self.map.insert(kmer, evidence);
    }

    /// Return `true` if the index contains any evidence for `kmer`.
    ///
    /// # Example
    /// ```
    /// use kam_index::HashKmerIndex;
    /// use kam_core::kmer::{KmerIndex, MoleculeEvidence};
    ///
    /// let mut index = HashKmerIndex::new();
    /// assert!(!index.contains(10));
    /// index.insert(10, MoleculeEvidence::default());
    /// assert!(index.contains(10));
    /// ```
    fn contains(&self, kmer: u64) -> bool {
        self.map.contains_key(&kmer)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use kam_core::kmer::KmerIndex;

    fn evidence(n_molecules: u32) -> MoleculeEvidence {
        MoleculeEvidence {
            n_molecules,
            ..Default::default()
        }
    }

    // Test 1: new index is empty
    #[test]
    fn new_index_is_empty() {
        let index = HashKmerIndex::new();
        assert!(index.is_empty());
        assert_eq!(index.len(), 0);
    }

    // Test 2: insert and get back same evidence
    #[test]
    fn insert_and_get_same_evidence() {
        let mut index = HashKmerIndex::new();
        let ev = MoleculeEvidence {
            n_molecules: 3,
            n_duplex: 2,
            n_simplex_fwd: 1,
            n_simplex_rev: 0,
            min_base_error_prob: 0.01,
            mean_base_error_prob: 0.02,
        };
        index.insert(42, ev.clone());
        let got = index.get(42).expect("evidence should be present");
        assert_eq!(*got, ev);
    }

    // Test 3: contains returns true after insert, false before
    #[test]
    fn contains_before_and_after_insert() {
        let mut index = HashKmerIndex::new();
        assert!(!index.contains(7));
        index.insert(7, evidence(1));
        assert!(index.contains(7));
        assert!(!index.contains(8));
    }

    // Test 4: molecule_count returns n_molecules from evidence
    #[test]
    fn molecule_count_returns_n_molecules() {
        let mut index = HashKmerIndex::new();
        index.insert(100, evidence(5));
        assert_eq!(index.molecule_count(100), 5);
    }

    // Test 5: molecule_count returns 0 for absent k-mer
    #[test]
    fn molecule_count_absent_returns_zero() {
        let index = HashKmerIndex::new();
        assert_eq!(index.molecule_count(999), 0);
    }

    // Test 6: entry() creates default if absent, returns existing if present
    #[test]
    fn entry_creates_default_and_returns_existing() {
        let mut index = HashKmerIndex::new();

        // First call creates a default entry
        {
            let ev = index.entry(55);
            assert_eq!(ev.n_molecules, 0);
            ev.n_molecules = 7;
        }

        // Second call returns the existing (modified) entry
        {
            let ev = index.entry(55);
            assert_eq!(ev.n_molecules, 7);
        }
    }

    // Test 7: multiple inserts to same k-mer: last write wins
    #[test]
    fn multiple_inserts_last_write_wins() {
        let mut index = HashKmerIndex::new();
        index.insert(1, evidence(10));
        index.insert(1, evidence(20));
        assert_eq!(index.molecule_count(1), 20);
        assert_eq!(index.len(), 1);
    }

    // Test 8: len() and is_empty() correct
    #[test]
    fn len_and_is_empty() {
        let mut index = HashKmerIndex::new();
        assert!(index.is_empty());
        assert_eq!(index.len(), 0);

        index.insert(1, evidence(1));
        assert!(!index.is_empty());
        assert_eq!(index.len(), 1);

        index.insert(2, evidence(2));
        assert_eq!(index.len(), 2);

        // Re-insert same key does not increase len
        index.insert(1, evidence(99));
        assert_eq!(index.len(), 2);
    }

    // Test 9: iter() yields all entries
    #[test]
    fn iter_yields_all_entries() {
        let mut index = HashKmerIndex::new();
        index.insert(10, evidence(1));
        index.insert(20, evidence(2));
        index.insert(30, evidence(3));

        let mut keys: Vec<u64> = index.iter().map(|(&k, _)| k).collect();
        keys.sort_unstable();
        assert_eq!(keys, vec![10, 20, 30]);

        let total_molecules: u32 = index.iter().map(|(_, ev)| ev.n_molecules).sum();
        assert_eq!(total_molecules, 6);
    }

    // Test 10: with_capacity doesn't change behaviour, just pre-allocates
    #[test]
    fn with_capacity_same_behaviour_as_new() {
        let mut index = HashKmerIndex::with_capacity(1_000);
        assert!(index.is_empty());
        assert_eq!(index.len(), 0);

        index.insert(42, evidence(3));
        assert!(index.contains(42));
        assert_eq!(index.molecule_count(42), 3);
        assert_eq!(index.len(), 1);
    }
}
