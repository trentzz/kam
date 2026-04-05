//! K-mer indexing types for the kam pipeline.
//!
//! The key insight is that [`MoleculeEvidence`] stores **molecule counts**, not
//! read counts. A k-mer present in 100 PCR duplicates of one molecule still has
//! `n_molecules = 1`. This is the fundamental advantage over Jellyfish.

/// Evidence associated with a single k-mer in the index.
///
/// All counts are at the **molecule** level. Two reads from the same molecule
/// both carrying a given k-mer contribute exactly `1` to `n_molecules`, not `2`.
///
/// # Example
/// ```
/// use kam_core::kmer::MoleculeEvidence;
///
/// let ev = MoleculeEvidence::default();
/// assert_eq!(ev.n_molecules, 0);
/// assert_eq!(ev.n_duplex, 0);
/// assert_eq!(ev.min_base_error_prob, 0.0);
/// ```
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct MoleculeEvidence {
    /// Number of distinct molecules that carry this k-mer.
    pub n_molecules: u32,
    /// Of `n_molecules`, how many have duplex support (both strands present).
    pub n_duplex: u32,
    /// Molecules with only forward-strand support.
    pub n_simplex_fwd: u32,
    /// Molecules with only reverse-strand support.
    pub n_simplex_rev: u32,
    /// Lowest per-base error probability observed across all supporting molecules
    /// (best quality observation; lower is better).
    pub min_base_error_prob: f32,
    /// Mean per-base error probability across all supporting molecules.
    pub mean_base_error_prob: f32,
}

/// An index mapping k-mers (encoded as `u64`) to [`MoleculeEvidence`].
///
/// Implementors include the in-memory hash map used by `kam-kmer`, as well as
/// any future on-disk or approximate (Bloom-filter) variants.
pub trait KmerIndex {
    /// Return a reference to the evidence for `kmer`, or `None` if absent.
    fn get(&self, kmer: u64) -> Option<&MoleculeEvidence>;

    /// Insert or overwrite evidence for `kmer`.
    fn insert(&mut self, kmer: u64, evidence: MoleculeEvidence);

    /// Return `true` if the index contains any evidence for `kmer`.
    fn contains(&self, kmer: u64) -> bool;

    /// Convenience method: return `n_molecules` for `kmer`, or `0` if absent.
    fn molecule_count(&self, kmer: u64) -> u32 {
        self.get(kmer).map_or(0, |ev| ev.n_molecules)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn molecule_evidence_default_all_zero() {
        let ev = MoleculeEvidence::default();
        assert_eq!(ev.n_molecules, 0);
        assert_eq!(ev.n_duplex, 0);
        assert_eq!(ev.n_simplex_fwd, 0);
        assert_eq!(ev.n_simplex_rev, 0);
        assert_eq!(ev.min_base_error_prob, 0.0);
        assert_eq!(ev.mean_base_error_prob, 0.0);
    }

    // A minimal in-memory implementation used only in tests.
    struct SimpleIndex(std::collections::HashMap<u64, MoleculeEvidence>);

    impl KmerIndex for SimpleIndex {
        fn get(&self, kmer: u64) -> Option<&MoleculeEvidence> {
            self.0.get(&kmer)
        }

        fn insert(&mut self, kmer: u64, evidence: MoleculeEvidence) {
            self.0.insert(kmer, evidence);
        }

        fn contains(&self, kmer: u64) -> bool {
            self.0.contains_key(&kmer)
        }
    }

    #[test]
    fn kmer_index_insert_get_contains() {
        let mut idx = SimpleIndex(std::collections::HashMap::new());
        let ev = MoleculeEvidence {
            n_molecules: 3,
            n_duplex: 2,
            ..Default::default()
        };
        idx.insert(42_u64, ev.clone());

        assert!(idx.contains(42));
        assert_eq!(idx.get(42).unwrap().n_molecules, 3);
        assert!(!idx.contains(99));
    }

    #[test]
    fn kmer_index_molecule_count_absent_returns_zero() {
        let idx = SimpleIndex(std::collections::HashMap::new());
        assert_eq!(idx.molecule_count(7), 0);
    }

    #[test]
    fn kmer_index_molecule_count_present() {
        let mut idx = SimpleIndex(std::collections::HashMap::new());
        idx.insert(
            1,
            MoleculeEvidence {
                n_molecules: 5,
                ..Default::default()
            },
        );
        assert_eq!(idx.molecule_count(1), 5);
    }
}
