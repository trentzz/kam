//! Fusion and translocation detection types and logic.
//!
//! A "fusion target" is a synthetic FASTA sequence formed by concatenating
//! breakpoint-adjacent segments from two partner loci. The FASTA header encodes
//! both partner coordinates so that a single entry carries all information needed
//! to call a BND record.
//!
//! The detection strategy reuses the existing walk/score/call pipeline without
//! modification. Fusion targets are walked like normal targets; the caller
//! handles the inverted ref/alt semantics (the fusion path IS the "reference"
//! for a fusion target — wild-type molecules do not span both loci).

use std::path::Path;

use kam_pathfind::score::PathEvidence;

use crate::caller::{assign_filter, compute_confidence, estimate_vaf, CallerConfig, VariantType};

// ─── Error type ───────────────────────────────────────────────────────────────

/// Errors that can occur when parsing fusion target definitions.
#[derive(Debug, thiserror::Error)]
pub enum FusionError {
    /// The FASTA header does not match the expected format.
    #[error("invalid fusion target header: {0}")]
    InvalidHeader(String),
    /// A genomic coordinate field could not be parsed.
    #[error("invalid coordinate in fusion header '{header}': {field}")]
    InvalidCoordinate {
        /// The full header string.
        header: String,
        /// Which field failed to parse.
        field: String,
    },
    /// The FASTA file could not be read.
    #[error("IO error reading fusion targets: {0}")]
    Io(#[from] std::io::Error),
    /// needletail returned a parse error.
    #[error("FASTA parse error: {0}")]
    Parse(String),
}

// ─── Public types ─────────────────────────────────────────────────────────────

/// A genomic coordinate range on a single chromosome.
///
/// # Example
/// ```
/// use kam_call::fusion::GenomicLocus;
/// let locus = GenomicLocus { chrom: "chr9".to_string(), start: 130_854_000, end: 130_854_050 };
/// assert_eq!(locus.chrom, "chr9");
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomicLocus {
    /// Chromosome name (e.g. "chr22").
    pub chrom: String,
    /// Start position (0-based, inclusive).
    pub start: u64,
    /// End position (0-based, exclusive).
    pub end: u64,
}

/// A synthetic fusion target representing the breakpoint junction between two loci.
///
/// The sequence is the concatenation of the partner A segment and the partner B
/// segment. The breakpoint is at `breakpoint_pos` (the number of bases from
/// partner A).
///
/// # Example
/// ```
/// use kam_call::fusion::{FusionTarget, GenomicLocus};
/// let target = FusionTarget {
///     name: "BCR_ABL1".to_string(),
///     locus_a: GenomicLocus { chrom: "chr22".to_string(), start: 23_632_500, end: 23_632_550 },
///     locus_b: GenomicLocus { chrom: "chr9".to_string(), start: 130_854_000, end: 130_854_050 },
///     sequence: b"ACGT".repeat(25).to_vec(),
///     breakpoint_pos: 50,
/// };
/// assert_eq!(target.name, "BCR_ABL1");
/// assert_eq!(target.breakpoint_pos, 50);
/// ```
#[derive(Debug, Clone)]
pub struct FusionTarget {
    /// Human-readable name (e.g. "BCR_ABL1").
    pub name: String,
    /// Genomic coordinates of the partner A (5') segment.
    pub locus_a: GenomicLocus,
    /// Genomic coordinates of the partner B (3') segment.
    pub locus_b: GenomicLocus,
    /// Concatenated sequence: partner A segment followed by partner B segment.
    pub sequence: Vec<u8>,
    /// Length of the partner A segment (index of the breakpoint in the sequence).
    pub breakpoint_pos: usize,
}

/// Context required to estimate VAF for a fusion call.
///
/// Wild-type molecules do not span both loci, so the fusion target's own depth
/// is not a valid denominator. Instead, the average depth at the two partner
/// loci provides the denominator.
///
/// # Example
/// ```
/// use kam_call::fusion::FusionContext;
/// let ctx = FusionContext { partner_a_depth: 1000.0, partner_b_depth: 900.0 };
/// assert_eq!(ctx.partner_a_depth, 1000.0);
/// ```
pub struct FusionContext {
    /// Mean molecule depth at partner A's normal target.
    pub partner_a_depth: f64,
    /// Mean molecule depth at partner B's normal target.
    pub partner_b_depth: f64,
}

/// A called fusion variant with full statistical evidence.
///
/// # Example
/// ```
/// use kam_call::fusion::{FusionCall, GenomicLocus};
/// use kam_call::caller::VariantFilter;
/// let call = FusionCall {
///     name: "BCR_ABL1".to_string(),
///     locus_a: GenomicLocus { chrom: "chr22".to_string(), start: 23_632_500, end: 23_632_550 },
///     locus_b: GenomicLocus { chrom: "chr9".to_string(), start: 130_854_000, end: 130_854_050 },
///     vaf: 0.01,
///     vaf_ci_low: 0.005,
///     vaf_ci_high: 0.02,
///     n_molecules: 10,
///     n_duplex: 5,
///     confidence: 0.999,
///     filter: VariantFilter::Pass,
/// };
/// assert_eq!(call.n_molecules, 10);
/// ```
#[derive(Debug, Clone)]
pub struct FusionCall {
    /// Human-readable name identifying this fusion (e.g. "BCR_ABL1").
    pub name: String,
    /// Genomic coordinates of the partner A breakpoint.
    pub locus_a: GenomicLocus,
    /// Genomic coordinates of the partner B breakpoint.
    pub locus_b: GenomicLocus,
    /// VAF point estimate.
    pub vaf: f64,
    /// Lower bound of the 95% Beta credible interval.
    pub vaf_ci_low: f64,
    /// Upper bound of the 95% Beta credible interval.
    pub vaf_ci_high: f64,
    /// Molecule count supporting the fusion.
    pub n_molecules: u32,
    /// Duplex molecule count supporting the fusion.
    pub n_duplex: u32,
    /// Posterior confidence that the fusion is real.
    pub confidence: f64,
    /// Quality filter outcome.
    pub filter: crate::caller::VariantFilter,
}

// ─── Public functions ─────────────────────────────────────────────────────────

/// Return true when `target_id` identifies a fusion target.
///
/// Fusion target IDs carry the suffix `__fusion` (two underscores, then
/// "fusion"). This sentinel distinguishes them from normal targets and SV
/// junction targets.
///
/// # Example
/// ```
/// use kam_call::fusion::is_fusion_target;
/// assert!(is_fusion_target("BCR_ABL1__chr22:23632500-23632550__chr9:130854000-130854050__fusion"));
/// assert!(!is_fusion_target("TP53_exon7"));
/// assert!(!is_fusion_target("chr5:1000-2000__del_junction"));
/// ```
pub fn is_fusion_target(target_id: &str) -> bool {
    target_id.ends_with("__fusion")
}

/// Parse a fusion target FASTA header into a [`FusionTarget`] (without sequence).
///
/// Expected format:
/// `{name}__{chromA}:{startA}-{endA}__{chromB}:{startB}-{endB}__fusion`
///
/// The sequence and `breakpoint_pos` are populated separately because the
/// header alone does not contain the sequence. `breakpoint_pos` is set to the
/// length of the partner A segment (`endA - startA`) as encoded in the header.
///
/// # Errors
///
/// Returns [`FusionError::InvalidHeader`] if the format does not match.
/// Returns [`FusionError::InvalidCoordinate`] if any numeric field is invalid.
///
/// # Example
/// ```
/// use kam_call::fusion::parse_fusion_header;
/// let target = parse_fusion_header(
///     "BCR_ABL1__chr22:23632500-23632550__chr9:130854000-130854050__fusion"
/// ).unwrap();
/// assert_eq!(target.name, "BCR_ABL1");
/// assert_eq!(target.locus_a.chrom, "chr22");
/// assert_eq!(target.locus_a.start, 23_632_500);
/// assert_eq!(target.locus_a.end, 23_632_550);
/// assert_eq!(target.locus_b.chrom, "chr9");
/// assert_eq!(target.locus_b.start, 130_854_000);
/// assert_eq!(target.locus_b.end, 130_854_050);
/// assert_eq!(target.breakpoint_pos, 50);
/// ```
pub fn parse_fusion_header(header: &str) -> Result<FusionTarget, FusionError> {
    // Strip the __fusion suffix.
    let body = header
        .strip_suffix("__fusion")
        .ok_or_else(|| FusionError::InvalidHeader(header.to_string()))?;

    // Split on double underscores: name, locusA, locusB.
    let parts: Vec<&str> = body.splitn(3, "__").collect();
    if parts.len() != 3 {
        return Err(FusionError::InvalidHeader(header.to_string()));
    }
    let name = parts[0].to_string();
    let locus_a = parse_locus(parts[1], header)?;
    let locus_b = parse_locus(parts[2], header)?;
    let breakpoint_pos = (locus_a.end - locus_a.start) as usize;

    Ok(FusionTarget {
        name,
        locus_a,
        locus_b,
        sequence: Vec::new(),
        breakpoint_pos,
    })
}

/// Load fusion targets from a FASTA file.
///
/// Each FASTA entry must follow the header format described in
/// [`parse_fusion_header`]. Only entries whose ID ends with `__fusion` are
/// loaded; other entries are silently skipped so that a combined targets file
/// (normal targets + fusion targets) is accepted without error.
///
/// # Errors
///
/// Returns [`FusionError::Io`] if the file cannot be opened.
/// Returns [`FusionError::Parse`] if needletail cannot parse the FASTA.
/// Returns [`FusionError::InvalidHeader`] or [`FusionError::InvalidCoordinate`]
/// if a fusion entry header is malformed.
///
/// # Example
/// ```no_run
/// use kam_call::fusion::parse_fusion_targets;
/// let targets = parse_fusion_targets(std::path::Path::new("fusion_targets.fa")).unwrap();
/// println!("loaded {} fusion targets", targets.len());
/// ```
pub fn parse_fusion_targets(path: &Path) -> Result<Vec<FusionTarget>, FusionError> {
    let mut reader =
        needletail::parse_fastx_file(path).map_err(|e| FusionError::Parse(e.to_string()))?;

    let mut targets = Vec::new();

    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| FusionError::Parse(e.to_string()))?;
        let id = std::str::from_utf8(rec.id())
            .map_err(|e| FusionError::Parse(e.to_string()))?
            .trim()
            .to_string();

        if !is_fusion_target(&id) {
            continue;
        }

        let mut target = parse_fusion_header(&id)?;
        target.sequence = rec.seq().to_vec();
        // Update breakpoint_pos based on actual sequence if locus coords match.
        // The header-derived breakpoint_pos (end_a - start_a) is the canonical
        // value; the actual sequence length of partner A is the same when the
        // FASTA was generated correctly.
        targets.push(target);
    }

    Ok(targets)
}

/// Call a fusion variant from the evidence on the fusion path.
///
/// The reference path through a fusion target IS the fusion allele. Wild-type
/// molecules do not produce k-mers on the fusion target, so the VAF denominator
/// must come from the partner loci depths rather than from the fusion target
/// itself.
///
/// When the partner depth denominator is zero (both partners have no coverage),
/// `None` is returned — there is no sensible call.
///
/// `ref_evidence` is the [`PathEvidence`] of the path that matched the fusion
/// target sequence (the fusion path, flagged as `is_reference` by
/// `score_and_rank_paths`).
///
/// # Example
/// ```
/// use kam_call::fusion::{FusionContext, FusionTarget, GenomicLocus, call_fusion};
/// use kam_call::caller::CallerConfig;
/// use kam_pathfind::score::PathEvidence;
///
/// let target = FusionTarget {
///     name: "BCR_ABL1".to_string(),
///     locus_a: GenomicLocus { chrom: "chr22".to_string(), start: 23_632_500, end: 23_632_550 },
///     locus_b: GenomicLocus { chrom: "chr9".to_string(), start: 130_854_000, end: 130_854_050 },
///     sequence: b"ACGT".repeat(25).to_vec(),
///     breakpoint_pos: 50,
/// };
/// let context = FusionContext { partner_a_depth: 1000.0, partner_b_depth: 900.0 };
/// let ev = PathEvidence {
///     min_molecules: 10,
///     mean_molecules: 10.0,
///     min_duplex: 5,
///     mean_duplex: 5.0,
///     min_variant_specific_duplex: 5,
///     mean_variant_specific_molecules: 10.0,
///     min_simplex_fwd: 5,
///     min_simplex_rev: 5,
///     mean_error_prob: 0.001,
/// };
/// let call = call_fusion(&ev, &context, &target, &CallerConfig::default());
/// assert!(call.is_some());
/// assert!(call.unwrap().n_molecules > 0);
/// ```
pub fn call_fusion(
    ref_evidence: &PathEvidence,
    context: &FusionContext,
    target: &FusionTarget,
    config: &CallerConfig,
) -> Option<FusionCall> {
    // Mean depth across both partners is the denominator for VAF.
    let partner_depth = (context.partner_a_depth + context.partner_b_depth) / 2.0;
    if partner_depth <= 0.0 {
        return None;
    }

    // Fusion molecule count: mean_variant_specific_molecules from the fusion path.
    // This mirrors the SV calling logic in call_variant.
    let k = ref_evidence.mean_variant_specific_molecules.round() as u32;
    let total = (partner_depth.round() as u32).saturating_add(k);

    let (vaf, vaf_ci_low, vaf_ci_high) = estimate_vaf(k, total);
    let confidence = compute_confidence(k, total, config.background_error_rate);

    let n_duplex = ref_evidence.min_variant_specific_duplex;

    let filter = assign_filter(
        confidence,
        1.0, // strand bias disabled for fusions (structurally asymmetric)
        k,
        n_duplex,
        vaf,
        VariantType::Fusion,
        config,
    );

    Some(FusionCall {
        name: target.name.clone(),
        locus_a: target.locus_a.clone(),
        locus_b: target.locus_b.clone(),
        vaf,
        vaf_ci_low,
        vaf_ci_high,
        n_molecules: k,
        n_duplex,
        confidence,
        filter,
    })
}

// ─── Private helpers ──────────────────────────────────────────────────────────

/// Parse a single genomic locus from the format `chrom:start-end`.
fn parse_locus(s: &str, header: &str) -> Result<GenomicLocus, FusionError> {
    // Split on ':' to separate chrom from the range.
    let colon = s
        .find(':')
        .ok_or_else(|| FusionError::InvalidHeader(header.to_string()))?;
    let chrom = s[..colon].to_string();
    let range = &s[colon + 1..];

    // Split range on '-'.
    let dash = range
        .find('-')
        .ok_or_else(|| FusionError::InvalidHeader(header.to_string()))?;
    let start_str = &range[..dash];
    let end_str = &range[dash + 1..];

    let start = start_str
        .parse::<u64>()
        .map_err(|_| FusionError::InvalidCoordinate {
            header: header.to_string(),
            field: format!("start={start_str:?}"),
        })?;
    let end = end_str
        .parse::<u64>()
        .map_err(|_| FusionError::InvalidCoordinate {
            header: header.to_string(),
            field: format!("end={end_str:?}"),
        })?;

    Ok(GenomicLocus { chrom, start, end })
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::caller::{CallerConfig, VariantFilter};
    use kam_pathfind::score::PathEvidence;

    fn make_evidence(mean_specific_molecules: f32, min_duplex: u32) -> PathEvidence {
        PathEvidence {
            min_molecules: mean_specific_molecules.round() as u32,
            mean_molecules: mean_specific_molecules,
            min_duplex,
            mean_duplex: min_duplex as f32,
            min_variant_specific_duplex: min_duplex,
            mean_variant_specific_molecules: mean_specific_molecules,
            min_simplex_fwd: 0,
            min_simplex_rev: 0,
            mean_error_prob: 0.001,
        }
    }

    fn bcr_abl1_target() -> FusionTarget {
        FusionTarget {
            name: "BCR_ABL1".to_string(),
            locus_a: GenomicLocus {
                chrom: "chr22".to_string(),
                start: 23_632_500,
                end: 23_632_550,
            },
            locus_b: GenomicLocus {
                chrom: "chr9".to_string(),
                start: 130_854_000,
                end: 130_854_050,
            },
            sequence: b"ACGT".repeat(25).to_vec(),
            breakpoint_pos: 50,
        }
    }

    // ── is_fusion_target ──────────────────────────────────────────────────────

    /// Fusion target ID with __fusion suffix is recognised.
    #[test]
    fn is_fusion_target_returns_true_for_fusion_id() {
        assert!(is_fusion_target(
            "BCR_ABL1__chr22:23632500-23632550__chr9:130854000-130854050__fusion"
        ));
    }

    /// Normal target IDs are not fusion targets.
    #[test]
    fn is_fusion_target_returns_false_for_normal_id() {
        assert!(!is_fusion_target("TP53_exon7"));
        assert!(!is_fusion_target("chr5:1000-2000"));
    }

    /// Empty string is not a fusion target.
    #[test]
    fn is_fusion_target_returns_false_for_empty() {
        assert!(!is_fusion_target(""));
    }

    /// A string that ends with "fusion" but not "__fusion" is not a fusion target.
    #[test]
    fn is_fusion_target_requires_double_underscore() {
        assert!(!is_fusion_target("BCR_ABL1_fusion"));
        assert!(!is_fusion_target("some_fusion"));
    }

    // ── parse_fusion_header ───────────────────────────────────────────────────

    /// A well-formed header parses correctly.
    #[test]
    fn parse_fusion_header_valid() {
        let target = parse_fusion_header(
            "BCR_ABL1__chr22:23632500-23632550__chr9:130854000-130854050__fusion",
        )
        .expect("should parse");
        assert_eq!(target.name, "BCR_ABL1");
        assert_eq!(target.locus_a.chrom, "chr22");
        assert_eq!(target.locus_a.start, 23_632_500);
        assert_eq!(target.locus_a.end, 23_632_550);
        assert_eq!(target.locus_b.chrom, "chr9");
        assert_eq!(target.locus_b.start, 130_854_000);
        assert_eq!(target.locus_b.end, 130_854_050);
        // breakpoint_pos = end_a - start_a = 50
        assert_eq!(target.breakpoint_pos, 50);
    }

    /// A header without the __fusion suffix is rejected.
    #[test]
    fn parse_fusion_header_no_suffix_fails() {
        let err =
            parse_fusion_header("BCR_ABL1__chr22:23632500-23632550__chr9:130854000-130854050")
                .expect_err("should fail");
        assert!(
            matches!(err, FusionError::InvalidHeader(_)),
            "expected InvalidHeader, got {err:?}"
        );
    }

    /// A header with only one locus field is rejected.
    #[test]
    fn parse_fusion_header_missing_locus_fails() {
        let err = parse_fusion_header("BCR_ABL1__chr22:23632500-23632550__fusion")
            .expect_err("should fail");
        assert!(matches!(
            err,
            FusionError::InvalidHeader(_) | FusionError::InvalidCoordinate { .. }
        ));
    }

    /// A non-numeric coordinate is rejected with an InvalidCoordinate error.
    #[test]
    fn parse_fusion_header_bad_coordinate_fails() {
        let err =
            parse_fusion_header("BCR_ABL1__chr22:ABC-23632550__chr9:130854000-130854050__fusion")
                .expect_err("should fail");
        assert!(
            matches!(err, FusionError::InvalidCoordinate { .. }),
            "expected InvalidCoordinate, got {err:?}"
        );
    }

    /// Locus with equal start and end parses (zero-length segment, unusual but valid).
    #[test]
    fn parse_fusion_header_zero_length_locus() {
        let target =
            parse_fusion_header("FGFR3_TACC3__chr4:1803548-1803548__chr4:1739323-1739373__fusion")
                .expect("should parse zero-length locus A");
        assert_eq!(target.locus_a.start, 1_803_548);
        assert_eq!(target.locus_a.end, 1_803_548);
        assert_eq!(target.breakpoint_pos, 0);
    }

    // ── call_fusion ───────────────────────────────────────────────────────────

    /// A call with good evidence passes all filters.
    #[test]
    fn call_fusion_pass_with_good_evidence() {
        let ev = make_evidence(10.0, 5);
        let ctx = FusionContext {
            partner_a_depth: 1000.0,
            partner_b_depth: 900.0,
        };
        let target = bcr_abl1_target();
        let config = CallerConfig::default();
        let call = call_fusion(&ev, &ctx, &target, &config).expect("should produce a call");
        assert_eq!(call.filter, VariantFilter::Pass);
        assert_eq!(call.n_molecules, 10);
        assert_eq!(call.n_duplex, 5);
        assert!(call.vaf > 0.0);
        assert!(call.confidence > 0.95);
    }

    /// Zero partner depth returns None.
    #[test]
    fn call_fusion_returns_none_for_zero_partner_depth() {
        let ev = make_evidence(10.0, 5);
        let ctx = FusionContext {
            partner_a_depth: 0.0,
            partner_b_depth: 0.0,
        };
        let target = bcr_abl1_target();
        let config = CallerConfig::default();
        let result = call_fusion(&ev, &ctx, &target, &config);
        assert!(result.is_none());
    }

    /// A call with zero fusion molecules is LowConfidence.
    #[test]
    fn call_fusion_low_confidence_for_zero_molecules() {
        let ev = make_evidence(0.0, 0);
        let ctx = FusionContext {
            partner_a_depth: 1000.0,
            partner_b_depth: 1000.0,
        };
        let target = bcr_abl1_target();
        let config = CallerConfig::default();
        let call = call_fusion(&ev, &ctx, &target, &config).expect("should produce a call");
        assert_eq!(call.filter, VariantFilter::LowConfidence);
        assert_eq!(call.n_molecules, 0);
    }

    /// VAF is estimated using the partner depth as the denominator.
    #[test]
    fn call_fusion_vaf_uses_partner_depth() {
        // 10 fusion molecules against a combined pool of 1000 from each partner.
        // mean_partner = 1000, so VAF ≈ 10 / (1000 + 10) ≈ 0.0099.
        let ev = make_evidence(10.0, 3);
        let ctx = FusionContext {
            partner_a_depth: 1000.0,
            partner_b_depth: 1000.0,
        };
        let target = bcr_abl1_target();
        let config = CallerConfig::default();
        let call = call_fusion(&ev, &ctx, &target, &config).expect("call");
        // VAF should be small (approximately 1%).
        assert!(call.vaf < 0.02, "expected VAF < 2%, got {}", call.vaf);
        assert!(call.vaf > 0.005, "expected VAF > 0.5%, got {}", call.vaf);
    }

    /// FusionCall fields are populated from the target.
    #[test]
    fn call_fusion_populates_loci_from_target() {
        let ev = make_evidence(5.0, 2);
        let ctx = FusionContext {
            partner_a_depth: 500.0,
            partner_b_depth: 500.0,
        };
        let target = bcr_abl1_target();
        let config = CallerConfig::default();
        let call = call_fusion(&ev, &ctx, &target, &config).expect("call");
        assert_eq!(call.name, "BCR_ABL1");
        assert_eq!(call.locus_a.chrom, "chr22");
        assert_eq!(call.locus_b.chrom, "chr9");
    }

    // ── parse_fusion_targets ──────────────────────────────────────────────────

    /// parse_fusion_targets skips non-fusion entries and loads fusion entries.
    #[test]
    fn parse_fusion_targets_from_fasta() {
        use std::io::Write;

        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("fusion.fa");
        {
            let mut f = std::fs::File::create(&path).expect("create");
            writeln!(
                f,
                ">BCR_ABL1__chr22:23632500-23632550__chr9:130854000-130854050__fusion"
            )
            .unwrap();
            writeln!(f, "{}", "A".repeat(100)).unwrap();
            // A normal (non-fusion) entry that should be skipped.
            writeln!(f, ">TP53_exon7").unwrap();
            writeln!(f, "{}", "C".repeat(100)).unwrap();
        }

        let targets = parse_fusion_targets(&path).expect("parse");
        assert_eq!(targets.len(), 1, "only fusion entries should be loaded");
        let t = &targets[0];
        assert_eq!(t.name, "BCR_ABL1");
        assert_eq!(t.sequence.len(), 100);
        assert_eq!(t.breakpoint_pos, 50);
    }

    /// parse_fusion_targets on a file with no fusion entries returns an empty vec.
    #[test]
    fn parse_fusion_targets_empty_when_no_fusion_entries() {
        use std::io::Write;

        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("targets.fa");
        {
            let mut f = std::fs::File::create(&path).expect("create");
            writeln!(f, ">TP53_exon7").unwrap();
            writeln!(f, "{}", "A".repeat(100)).unwrap();
        }

        let targets = parse_fusion_targets(&path).expect("parse");
        assert!(targets.is_empty());
    }

    // ── Additional edge-case tests ───────────────────────────────────────────

    // Test: parse_fusion_header with completely malformed header (no double underscores).
    // Must return an error, not panic.
    #[test]
    fn parse_fusion_header_completely_malformed() {
        let err = parse_fusion_header("just_a_random_string").expect_err("should fail");
        assert!(
            matches!(err, FusionError::InvalidHeader(_)),
            "expected InvalidHeader, got {err:?}"
        );
    }

    // Test: parse_fusion_header with missing end coordinate.
    // "chr22:23632500-" has a dash but no end value, which should fail to parse.
    #[test]
    fn parse_fusion_header_missing_end_coordinate() {
        let err = parse_fusion_header(
            "BCR_ABL1__chr22:23632500-__chr9:130854000-130854050__fusion",
        )
        .expect_err("should fail");
        assert!(
            matches!(err, FusionError::InvalidCoordinate { .. }),
            "expected InvalidCoordinate, got {err:?}"
        );
    }

    // Test: parse_fusion_targets on an empty FASTA file returns an empty vec without error.
    // An empty file is a valid edge case in production (no targets configured).
    #[test]
    fn parse_fusion_targets_empty_file() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("empty.fa");
        {
            let _f = std::fs::File::create(&path).expect("create");
            // Write nothing; the file is empty.
        }
        // needletail on an empty file may return an error or empty; either is acceptable.
        // The important thing is it does not panic.
        let result = parse_fusion_targets(&path);
        if let Ok(targets) = result {
            assert!(targets.is_empty(), "empty file should produce no targets");
        }
        // An IO or parse error on an empty file is also acceptable.
    }

    // Test: call_fusion with zero partner depth on one side and positive on the other.
    // The mean partner depth is positive, so a call should still be produced.
    #[test]
    fn call_fusion_one_partner_zero_depth() {
        let ev = make_evidence(5.0, 2);
        let ctx = FusionContext {
            partner_a_depth: 0.0,
            partner_b_depth: 1000.0,
        };
        let target = bcr_abl1_target();
        let config = CallerConfig::default();
        let call = call_fusion(&ev, &ctx, &target, &config);
        // partner_depth = (0 + 1000) / 2 = 500, which is > 0, so a call is produced.
        assert!(
            call.is_some(),
            "one zero partner should still allow a call when the other is positive"
        );
    }

    // Test: call_fusion confidence interval is ordered (low <= vaf <= high).
    // This is a basic invariant that must hold for all valid calls.
    #[test]
    fn call_fusion_ci_is_ordered() {
        let ev = make_evidence(8.0, 3);
        let ctx = FusionContext {
            partner_a_depth: 800.0,
            partner_b_depth: 800.0,
        };
        let target = bcr_abl1_target();
        let config = CallerConfig::default();
        let call = call_fusion(&ev, &ctx, &target, &config).expect("should produce a call");
        assert!(
            call.vaf_ci_low <= call.vaf,
            "CI lower ({}) must be <= VAF ({})",
            call.vaf_ci_low,
            call.vaf
        );
        assert!(
            call.vaf <= call.vaf_ci_high,
            "VAF ({}) must be <= CI upper ({})",
            call.vaf,
            call.vaf_ci_high
        );
    }
}
