//! Core types for variant calling: enums, structs, and their trait implementations.

// ─── CallSource ───────────────────────────────────────────────────────────────

/// Indicates how a variant call record was produced.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, serde::Serialize)]
pub enum CallSource {
    /// Normal pathfinding call that passed or failed statistical thresholds.
    #[default]
    Called,
    /// Pathfinding found a matching path but confidence was below threshold.
    SubThreshold,
    /// No call from pathfinding; k-mer index evidence found at target position.
    Rescued,
    /// No call from pathfinding; no k-mer evidence at target position.
    NoEvidence,
}

impl std::fmt::Display for CallSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CallSource::Called => write!(f, "CALLED"),
            CallSource::SubThreshold => write!(f, "SUBTHRESHOLD"),
            CallSource::Rescued => write!(f, "RESCUED"),
            CallSource::NoEvidence => write!(f, "NO_EVIDENCE"),
        }
    }
}

// ─── VariantType ──────────────────────────────────────────────────────────────

/// Minimum indel length to classify as a structural variant.
pub const SV_LENGTH_THRESHOLD: usize = 50;

/// Classification of a variant by allele length comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum VariantType {
    /// Single nucleotide variant.
    Snv,
    /// Insertion (alt longer than ref, < 50 bp).
    Insertion,
    /// Deletion (alt shorter than ref, < 50 bp).
    Deletion,
    /// Multi-nucleotide variant (same length, >1 difference).
    Mnv,
    /// Complex rearrangement.
    Complex,
    /// Large deletion (alt shorter than ref by ≥ 50 bp).
    LargeDeletion,
    /// Tandem duplication (alt longer than ref by ≥ 50 bp).
    TandemDuplication,
    /// Inversion (same length, alt is reverse-complement of ref segment).
    Inversion,
    /// Fusion: junction between two distinct genomic loci.
    ///
    /// Detected when the alt path connects k-mers from non-contiguous
    /// reference regions, implying a chromosomal translocation or gene fusion.
    Fusion,
    /// Inversion with flanking deletion.
    ///
    /// The alt sequence length differs from ref AND contains a central segment
    /// that is the reverse complement of the corresponding ref region. Both
    /// the length change and the RC segment must meet the SV length threshold.
    InvDel,
    /// Novel insertion of ≥ 50 bp that is not a tandem repeat of nearby
    /// reference sequence.
    ///
    /// Distinguished from `TandemDuplication` by the absence of sequence
    /// similarity between the inserted bases and the flanking reference.
    NovelInsertion,
}

impl std::fmt::Display for VariantType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            VariantType::Snv => "SNV",
            VariantType::Insertion => "Insertion",
            VariantType::Deletion => "Deletion",
            VariantType::Mnv => "MNV",
            VariantType::Complex => "Complex",
            VariantType::LargeDeletion => "LargeDeletion",
            VariantType::TandemDuplication => "TandemDuplication",
            VariantType::Inversion => "Inversion",
            VariantType::Fusion => "Fusion",
            VariantType::InvDel => "InvDel",
            VariantType::NovelInsertion => "NovelInsertion",
        };
        f.write_str(s)
    }
}

// ─── VariantFilter ────────────────────────────────────────────────────────────

/// Filter outcome for a variant call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantFilter {
    /// Variant passes all filters.
    Pass,
    /// Significant strand bias detected.
    StrandBias,
    /// Posterior confidence below threshold.
    LowConfidence,
    /// Fewer duplex molecules than required at the variant site.
    LowDuplex,
    /// UMI collision risk too high to trust.
    CollisionRisk,
    /// VAF exceeds the configured maximum (likely germline).
    HighVaf,
    /// Variant does not match any entry in the target variants set.
    ///
    /// Applied in tumour-informed monitoring mode (`--target-variants`).
    /// A call at this locus passed all statistical filters but the
    /// called allele is not one of the expected somatic variants.
    NotTargeted,
}

// ─── VariantCall ─────────────────────────────────────────────────────────────

/// A called variant with full statistical evidence.
///
/// # Example
/// ```
/// use kam_call::caller::{VariantFilter, VariantType};
/// let vt = VariantType::Snv;
/// let vf = VariantFilter::Pass;
/// assert_eq!(vt, VariantType::Snv);
/// assert_eq!(vf, VariantFilter::Pass);
/// ```
#[derive(Debug, Clone)]
pub struct VariantCall {
    /// Identifier for the target region.
    pub target_id: String,
    /// Variant class inferred from ref/alt sequence comparison.
    pub variant_type: VariantType,
    /// Reference allele sequence.
    pub ref_sequence: Vec<u8>,
    /// Alternate allele sequence.
    pub alt_sequence: Vec<u8>,
    /// VAF point estimate: `k / M`.
    pub vaf: f64,
    /// Lower bound of 95% Beta credible interval.
    pub vaf_ci_low: f64,
    /// Upper bound of 95% Beta credible interval.
    pub vaf_ci_high: f64,
    /// Number of molecules supporting the reference allele.
    pub n_molecules_ref: u32,
    /// Number of molecules supporting the alternate allele.
    pub n_molecules_alt: u32,
    /// Number of duplex molecules supporting the alternate allele.
    pub n_duplex_alt: u32,
    /// Number of simplex molecules supporting the alternate allele.
    pub n_simplex_alt: u32,
    /// Forward-strand simplex molecules supporting the alternate allele.
    ///
    /// Derived from `alt PathEvidence.min_simplex_fwd`.
    pub n_simplex_fwd_alt: u32,
    /// Reverse-strand simplex molecules supporting the alternate allele.
    ///
    /// Derived from `alt PathEvidence.min_simplex_rev`.
    pub n_simplex_rev_alt: u32,
    /// Duplex molecules supporting the reference allele.
    ///
    /// Derived from `ref PathEvidence.min_duplex`.
    pub n_duplex_ref: u32,
    /// Simplex molecules supporting the reference allele (both strands combined).
    ///
    /// Sum of `ref PathEvidence.min_simplex_fwd` and `ref PathEvidence.min_simplex_rev`.
    pub n_simplex_ref: u32,
    /// Mean per-base error probability for the alternate path.
    ///
    /// Derived from `alt PathEvidence.mean_error_prob`.
    pub mean_alt_error_prob: f32,
    /// Minimum duplex support at variant-specific k-mers in the alternate path.
    ///
    /// Derived from `alt PathEvidence.min_variant_specific_duplex`.
    pub min_variant_specific_duplex: u32,
    /// Mean molecule count at variant-specific k-mers in the alternate path.
    ///
    /// Derived from `alt PathEvidence.mean_variant_specific_molecules`.
    pub mean_variant_specific_molecules: f32,
    /// Posterior probability that the variant is real.
    pub confidence: f64,
    /// Fisher's exact test p-value for strand bias.
    pub strand_bias_p: f64,
    /// Quality filter outcome.
    pub filter: VariantFilter,
    /// ML model posterior probability (class 1 = real variant).
    ///
    /// `None` when no model was loaded.
    pub ml_prob: Option<f32>,
    /// How this call record was produced.
    pub call_source: CallSource,
    /// Minimum alt-specific k-mer molecule count (rescue probe only).
    pub rescue_min_alt_molecules: Option<u32>,
    /// Duplex count across alt-specific k-mers (rescue probe only).
    pub rescue_alt_duplex: Option<u32>,
    /// Approximate VAF from k-mer evidence (rescue probe only).
    pub rescue_approx_vaf: Option<f32>,
    /// Number of alt-specific k-mers found in index (rescue probe only).
    pub rescue_kmers_found: Option<usize>,
    /// Total alt-specific k-mers checked (rescue probe only).
    pub rescue_kmers_total: Option<usize>,
}
