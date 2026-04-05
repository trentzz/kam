//! Core molecule types for the kam pipeline.
//!
//! The [`Molecule`] is the atomic unit of evidence. All counts, grouping, and
//! evidence are molecule-level, never read-level.

use crate::kmer::MoleculeEvidence;

/// A consensus sequence derived from a set of reads from the same strand of the
/// same DNA molecule.
///
/// `per_base_error_prob` is the probability of error at each position (not Phred
/// quality — a value of `0.001` means 0.1% error probability).
///
/// `per_base_strand_support` counts how many forward-strand and reverse-strand
/// reads agree with the called base at each position.
///
/// `family_size` is `(fwd_reads, rev_reads)` — the number of reads from each
/// strand that contributed to this consensus.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ConsensusRead {
    /// The consensus base sequence (ASCII bytes, e.g. `b"ACGT..."`).
    pub sequence: Vec<u8>,
    /// Per-base probability of error (0.0–1.0).
    pub per_base_error_prob: Vec<f32>,
    /// Per-base strand support as `(fwd_count, rev_count)`.
    pub per_base_strand_support: Vec<(u8, u8)>,
    /// `(fwd_reads, rev_reads)` — total reads on each strand.
    pub family_size: (u8, u8),
}

/// A canonical, strand-agnostic identifier for a duplex UMI pair.
///
/// Two reads from opposite strands of the same molecule will have their UMIs
/// swapped relative to each other. `CanonicalUmiPair` resolves this by always
/// placing the lexicographically smaller UMI in `umi_a`, so both strands map to
/// the identical key.
///
/// UMI length is variable: any length is accepted. The standard Twist UMI
/// duplex chemistry uses 5 bp UMIs, but longer or shorter UMIs are supported.
///
/// # Example
/// ```
/// use kam_core::molecule::{CanonicalUmiPair, Strand};
///
/// // "ACGTA" < "TGCAT", so umi_a stays as-is.
/// let pair = CanonicalUmiPair::new(b"ACGTA".to_vec(), b"TGCAT".to_vec());
/// assert_eq!(pair.umi_a, b"ACGTA");
/// assert_eq!(pair.umi_b, b"TGCAT");
///
/// // Swapped inputs produce the identical canonical pair.
/// let pair2 = CanonicalUmiPair::new(b"TGCAT".to_vec(), b"ACGTA".to_vec());
/// assert_eq!(pair, pair2);
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub struct CanonicalUmiPair {
    /// The lexicographically smaller UMI.
    pub umi_a: Vec<u8>,
    /// The lexicographically larger UMI (or equal, for self-complementary UMIs).
    pub umi_b: Vec<u8>,
}

impl CanonicalUmiPair {
    /// Create a canonical UMI pair from R1 and R2 UMIs.
    ///
    /// The lexicographically smaller UMI is always stored as `umi_a`.
    /// This ensures forward and reverse strand reads from the same molecule
    /// hash to the same key.
    ///
    /// Accepts any UMI length; the standard Twist UMI duplex chemistry uses
    /// 5 bp UMIs.
    ///
    /// # Example
    /// ```
    /// use kam_core::molecule::CanonicalUmiPair;
    ///
    /// let pair = CanonicalUmiPair::new(b"ACGTA".to_vec(), b"TGCAT".to_vec());
    /// assert_eq!(pair.umi_a, b"ACGTA");
    /// assert_eq!(pair.umi_b, b"TGCAT");
    /// ```
    pub fn new(umi_r1: Vec<u8>, umi_r2: Vec<u8>) -> Self {
        if umi_r1 <= umi_r2 {
            Self {
                umi_a: umi_r1,
                umi_b: umi_r2,
            }
        } else {
            Self {
                umi_a: umi_r2,
                umi_b: umi_r1,
            }
        }
    }

    /// Determine which strand R1 came from based on whether its UMI is the
    /// canonical `umi_a` ([`Strand::Forward`]) or `umi_b` ([`Strand::Reverse`]).
    ///
    /// # Example
    /// ```
    /// use kam_core::molecule::{CanonicalUmiPair, Strand};
    ///
    /// let pair = CanonicalUmiPair::new(b"ACGTA".to_vec(), b"TGCAT".to_vec());
    /// assert_eq!(pair.strand_of_r1(b"ACGTA"), Strand::Forward);
    /// assert_eq!(pair.strand_of_r1(b"TGCAT"), Strand::Reverse);
    /// ```
    pub fn strand_of_r1(&self, umi_r1: &[u8]) -> Strand {
        if umi_r1 == self.umi_a {
            Strand::Forward
        } else {
            Strand::Reverse
        }
    }
}

/// Which strand of the original double-stranded DNA molecule a read came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Strand {
    /// The read came from the forward strand (R1 UMI is the canonical `umi_a`).
    Forward,
    /// The read came from the reverse strand (R1 UMI is the canonical `umi_b`).
    Reverse,
}

/// The type of read family, determined by how many reads are present on each
/// strand.
///
/// Duplex support (reads on both strands) is the strongest evidence class;
/// simplex and singleton reads are progressively weaker.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum FamilyType {
    /// At least one forward-strand read **and** at least one reverse-strand
    /// read are present.
    Duplex,
    /// Only forward-strand reads are present (no reverse-strand reads).
    SimplexFwd,
    /// Only reverse-strand reads are present (no forward-strand reads).
    SimplexRev,
    /// Exactly one read total (degenerate family).
    Singleton,
}

impl FamilyType {
    /// Classify a family given its `(fwd_reads, rev_reads)` family size.
    ///
    /// A `(0, 0)` family (no reads on either strand) returns `Singleton` as the
    /// most conservative fallback.
    ///
    /// # Example
    /// ```
    /// use kam_core::molecule::FamilyType;
    ///
    /// assert_eq!(FamilyType::from_family_size((2, 3)), FamilyType::Duplex);
    /// assert_eq!(FamilyType::from_family_size((1, 0)), FamilyType::Singleton);
    /// assert_eq!(FamilyType::from_family_size((3, 0)), FamilyType::SimplexFwd);
    /// assert_eq!(FamilyType::from_family_size((0, 2)), FamilyType::SimplexRev);
    /// assert_eq!(FamilyType::from_family_size((0, 0)), FamilyType::Singleton);
    /// ```
    pub fn from_family_size(family_size: (u8, u8)) -> Self {
        let (fwd, rev) = family_size;
        match (fwd, rev) {
            (0, 0) => FamilyType::Singleton,
            (1, 0) | (0, 1) => FamilyType::Singleton,
            (f, 0) if f > 0 => FamilyType::SimplexFwd,
            (0, r) if r > 0 => FamilyType::SimplexRev,
            _ => FamilyType::Duplex,
        }
    }
}

/// A fully-assembled molecule carrying consensus evidence from both strands.
///
/// Fields are populated by `kam-assemble` during UMI grouping and consensus
/// calling. Matches the `core_data_model` spec.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Molecule {
    /// Stable molecule identifier derived from hashing the canonical UMI pair.
    pub id: u64,
    /// Forward-strand UMI (variable length; 5 bp for Twist UMI duplex chemistry).
    pub umi_fwd: Vec<u8>,
    /// Reverse-strand UMI (variable length; 5 bp for Twist UMI duplex chemistry).
    pub umi_rev: Vec<u8>,
    /// Consensus from the forward-strand reads, if available.
    pub consensus_fwd: Option<ConsensusRead>,
    /// Consensus from the reverse-strand reads, if available.
    pub consensus_rev: Option<ConsensusRead>,
    /// Duplex consensus (both strands must agree), if available.
    pub duplex_consensus: Option<ConsensusRead>,
    /// Aggregate k-mer evidence for this molecule (populated by `kam-kmer`).
    pub evidence: Option<MoleculeEvidence>,
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── CanonicalUmiPair::new() ──────────────────────────────────────────────

    #[test]
    fn canonical_pair_r1_smaller() {
        let pair = CanonicalUmiPair::new(b"ACGTA".to_vec(), b"TGCAT".to_vec());
        assert_eq!(pair.umi_a, b"ACGTA");
        assert_eq!(pair.umi_b, b"TGCAT");
    }

    #[test]
    fn canonical_pair_r1_larger_swapped() {
        let pair = CanonicalUmiPair::new(b"TGCAT".to_vec(), b"ACGTA".to_vec());
        assert_eq!(pair.umi_a, b"ACGTA");
        assert_eq!(pair.umi_b, b"TGCAT");
    }

    #[test]
    fn canonical_pair_equal_umis() {
        let pair = CanonicalUmiPair::new(b"AAAAA".to_vec(), b"AAAAA".to_vec());
        assert_eq!(pair.umi_a, b"AAAAA");
        assert_eq!(pair.umi_b, b"AAAAA");
    }

    #[test]
    fn canonical_pair_all_a_vs_all_t() {
        // "AAAAA" < "TTTTT" lexicographically
        let pair = CanonicalUmiPair::new(b"TTTTT".to_vec(), b"AAAAA".to_vec());
        assert_eq!(pair.umi_a, b"AAAAA");
        assert_eq!(pair.umi_b, b"TTTTT");
    }

    #[test]
    fn canonical_pair_variable_length() {
        // 8 bp UMIs — not tied to the Twist 5 bp preset.
        let pair = CanonicalUmiPair::new(b"ACGTACGT".to_vec(), b"TGCATGCA".to_vec());
        assert_eq!(pair.umi_a, b"ACGTACGT");
        assert_eq!(pair.umi_b, b"TGCATGCA");
    }

    // ── strand_of_r1() ──────────────────────────────────────────────────────

    #[test]
    fn strand_of_r1_forward() {
        let pair = CanonicalUmiPair::new(b"ACGTA".to_vec(), b"TGCAT".to_vec());
        assert_eq!(pair.strand_of_r1(b"ACGTA"), Strand::Forward);
    }

    #[test]
    fn strand_of_r1_reverse() {
        let pair = CanonicalUmiPair::new(b"ACGTA".to_vec(), b"TGCAT".to_vec());
        assert_eq!(pair.strand_of_r1(b"TGCAT"), Strand::Reverse);
    }

    #[test]
    fn same_molecule_both_orientations_produce_identical_pair() {
        // Read from forward strand: R1=ACGTA, R2=TGCAT
        let pair_fwd = CanonicalUmiPair::new(b"ACGTA".to_vec(), b"TGCAT".to_vec());
        // Read from reverse strand: R1 and R2 are swapped
        let pair_rev = CanonicalUmiPair::new(b"TGCAT".to_vec(), b"ACGTA".to_vec());
        assert_eq!(pair_fwd, pair_rev);
    }

    #[test]
    fn canonical_pair_hash_consistent() {
        use std::collections::HashSet;
        let pair_a = CanonicalUmiPair::new(b"ACGTA".to_vec(), b"TGCAT".to_vec());
        let pair_b = CanonicalUmiPair::new(b"TGCAT".to_vec(), b"ACGTA".to_vec());
        let mut set = HashSet::new();
        set.insert(pair_a.clone());
        // The same canonical pair (built from swapped inputs) must already be in the set.
        assert!(set.contains(&pair_b));
    }

    // ── FamilyType ───────────────────────────────────────────────────────────

    #[test]
    fn family_type_duplex() {
        assert_eq!(FamilyType::from_family_size((2, 3)), FamilyType::Duplex);
    }

    #[test]
    fn family_type_singleton_fwd() {
        assert_eq!(FamilyType::from_family_size((1, 0)), FamilyType::Singleton);
    }

    #[test]
    fn family_type_singleton_rev() {
        assert_eq!(FamilyType::from_family_size((0, 1)), FamilyType::Singleton);
    }

    #[test]
    fn family_type_simplex_fwd() {
        assert_eq!(FamilyType::from_family_size((3, 0)), FamilyType::SimplexFwd);
    }

    #[test]
    fn family_type_simplex_rev() {
        assert_eq!(FamilyType::from_family_size((0, 2)), FamilyType::SimplexRev);
    }

    #[test]
    fn family_type_zero_zero_is_singleton() {
        // A family with no reads on either strand must not be classified as
        // Duplex. The catch-all arm previously matched this case incorrectly.
        assert_eq!(FamilyType::from_family_size((0, 0)), FamilyType::Singleton);
    }
}
