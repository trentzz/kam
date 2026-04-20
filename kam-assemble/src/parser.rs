//! UMI FASTQ read pair parser.
//!
//! Extracts UMI, skip bases, and template from raw R1/R2 reads according to
//! the configured [`ReadStructure`].  The parser is chemistry-agnostic: it
//! accepts any UMI length and any skip length (including 0 bp).  It is also
//! k-agnostic: it does **not** enforce a minimum template length by default.
//! Short templates are parsed successfully; downstream stages (e.g. k-mer
//! extraction) decide whether to use them.

use kam_core::chemistry::ReadStructure;
use kam_core::molecule::{CanonicalUmiPair, Strand};

// ── Configuration ────────────────────────────────────────────────────────────

/// Configuration passed to [`parse_read_pair`].
///
/// Use [`ParserConfig::default`] for the standard Twist UMI duplex preset with
/// no quality or length filters, or set the optional fields to add filters.
#[derive(Debug, Clone)]
pub struct ParserConfig {
    /// Describes where in each read the UMI, skip, and template regions live.
    pub read_structure: ReadStructure,
    /// Reject templates shorter than this many bases.
    ///
    /// `None` (the default) accepts templates of any length, including 0 bp.
    /// This keeps the parser k-agnostic — the k-mer extraction stage decides
    /// the effective minimum.
    pub min_template_length: Option<usize>,
    /// Reject UMIs that contain any base with Phred quality below this
    /// threshold.
    ///
    /// `None` (the default) accepts UMIs regardless of base quality.
    pub min_umi_quality: Option<u8>,
}

impl Default for ParserConfig {
    /// Default config: Twist UMI duplex chemistry, no filters.
    fn default() -> Self {
        Self {
            read_structure: ReadStructure::twist_umi_duplex(),
            min_template_length: None,
            min_umi_quality: None,
        }
    }
}

// ── Drop reason ───────────────────────────────────────────────────────────────

/// Why a read pair was rejected by [`parse_read_pair`].
///
/// Carried inside [`ParseResult::Dropped`] so that callers can increment the
/// right counter in [`ParseStats`] and optionally write a structured drop log.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DropReason {
    /// The read is shorter than `umi_length + skip_length` — not even enough
    /// bytes to extract the UMI and skip.  This usually indicates a truncated
    /// or malformed record.
    ReadTooShort,
    /// The template region is shorter than the user-configured
    /// [`ParserConfig::min_template_length`].
    TemplateTooShort,
    /// At least one UMI base has a Phred quality below
    /// [`ParserConfig::min_umi_quality`].
    LowUmiQuality,
}

// ── Output types ─────────────────────────────────────────────────────────────

/// All fields extracted from a successfully parsed read pair.
#[derive(Debug, Clone)]
pub struct ParsedReadPair {
    /// UMI bases from R1 (positions 0..umi_length).
    pub umi_r1: Vec<u8>,
    /// UMI bases from R2 (positions 0..umi_length).
    pub umi_r2: Vec<u8>,
    /// Skip bases from R1 (positions umi_length..template_start).
    pub skip_r1: Vec<u8>,
    /// Skip bases from R2 (positions umi_length..template_start).
    pub skip_r2: Vec<u8>,
    /// Template sequence from R1 (positions template_start..).
    pub template_r1: Vec<u8>,
    /// Template sequence from R2 (positions template_start..).
    pub template_r2: Vec<u8>,
    /// Base qualities corresponding to `template_r1`.
    pub qual_r1: Vec<u8>,
    /// Base qualities corresponding to `template_r2`.
    pub qual_r2: Vec<u8>,
    /// Base qualities for the UMI bases of R1.
    pub umi_qual_r1: Vec<u8>,
    /// Base qualities for the UMI bases of R2.
    pub umi_qual_r2: Vec<u8>,
    /// Strand-agnostic canonical UMI pair derived from `umi_r1` and `umi_r2`.
    pub canonical_umi: CanonicalUmiPair,
    /// Which strand R1 came from, as determined by the canonical pairing.
    pub strand: Strand,
}

/// The outcome of [`parse_read_pair`].
///
/// Both variants carry enough information for structured logging: `Ok` holds
/// all extracted fields, `Dropped` carries the machine-readable [`DropReason`]
/// and a human-readable `detail` string.
///
/// The `Ok` variant boxes its payload to keep the enum size small.
#[derive(Debug)]
pub enum ParseResult {
    /// The read pair was parsed successfully.
    Ok(Box<ParsedReadPair>),
    /// The read pair was rejected.
    Dropped {
        /// Machine-readable rejection category.
        reason: DropReason,
        /// Human-readable detail for log sinks (e.g. `"r1_len=4,min=7"`).
        detail: String,
    },
}

// ── Stats ─────────────────────────────────────────────────────────────────────

/// Running counters for a parsing session.
///
/// Intended to be maintained by the caller, not by [`parse_read_pair`] itself.
/// Counters are always incremented; detailed per-read log lines are only
/// emitted when the log sink is enabled (see `docs/planning/logging_architecture.md`).
///
/// # Example
/// ```
/// use kam_assemble::parser::ParseStats;
///
/// let mut stats = ParseStats::default();
/// stats.n_processed += 1;
/// stats.n_passed += 1;
/// ```
#[derive(Debug, Clone, Default)]
pub struct ParseStats {
    /// Total read pairs examined.
    pub n_processed: u64,
    /// Read pairs that produced a [`ParseResult::Ok`].
    pub n_passed: u64,
    /// Read pairs dropped because R1 or R2 was shorter than `umi+skip`.
    pub n_read_too_short: u64,
    /// Read pairs dropped because the template was shorter than
    /// [`ParserConfig::min_template_length`].
    pub n_template_too_short: u64,
    /// Read pairs dropped because a UMI base quality was below
    /// [`ParserConfig::min_umi_quality`].
    pub n_low_umi_quality: u64,
}

// ── Core parsing function ─────────────────────────────────────────────────────

/// Parse one R1/R2 read pair according to `config`.
///
/// Extracts the UMI (first `umi_length` bases), the skip region
/// (`umi_length..template_start`), and the template (`template_start..`)
/// from each read.  Quality arrays are split at the same positions.
///
/// Both UMI length and skip length are taken from [`ParserConfig::read_structure`]
/// at runtime.  Any skip length is accepted, including 0 bp (as used by the
/// `simplex_12bp` preset).
///
/// The parser applies filters in the following order:
/// 1. **ReadTooShort** — R1 or R2 is shorter than `umi_length + skip_length`.
/// 2. **TemplateTooShort** — template is shorter than `min_template_length`
///    (only checked when that option is `Some`).
/// 3. **LowUmiQuality** — any UMI base quality is below `min_umi_quality`
///    (only checked when that option is `Some`).
///
/// # Arguments
///
/// * `r1_seq` / `r1_qual` — raw base and quality bytes for R1.
/// * `r2_seq` / `r2_qual` — raw base and quality bytes for R2.
/// * `config` — parser configuration (chemistry preset + optional filters).
///
/// # Examples
///
/// ```
/// use kam_assemble::parser::{ParserConfig, ParseResult, parse_read_pair};
///
/// let config = ParserConfig::default();
/// // 5 bp UMI + 2 bp skip + 10 bp template = 17 bp total
/// let r1_seq  = b"ACGTATGNNNNNNNNNN";
/// let r1_qual = b"IIIIIIIIIIIIIIIII";
/// let r2_seq  = b"TGCATAGNNNNNNNNNN";
/// let r2_qual = b"IIIIIIIIIIIIIIIII";
/// let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &config);
/// assert!(matches!(result, ParseResult::Ok(_)));
/// ```
pub fn parse_read_pair(
    r1_seq: &[u8],
    r1_qual: &[u8],
    r2_seq: &[u8],
    r2_qual: &[u8],
    config: &ParserConfig,
) -> ParseResult {
    let umi_len = config.read_structure.umi_length;
    let tmpl_start = config.read_structure.template_start();

    // ── 1. ReadTooShort ───────────────────────────────────────────────────
    if r1_seq.len() < tmpl_start || r2_seq.len() < tmpl_start {
        let detail = format!(
            "r1_len={},r2_len={},min={}",
            r1_seq.len(),
            r2_seq.len(),
            tmpl_start
        );
        return ParseResult::Dropped {
            reason: DropReason::ReadTooShort,
            detail,
        };
    }

    // ── 2. Extract UMI, skip, template ───────────────────────────────────
    // Both UMI and skip use Vec<u8> so any length (including 0 bp skip) is
    // accepted without coercion or error.
    let umi_r1 = r1_seq[..umi_len].to_vec();
    let umi_r2 = r2_seq[..umi_len].to_vec();

    let skip_r1 = r1_seq[umi_len..tmpl_start].to_vec();
    let skip_r2 = r2_seq[umi_len..tmpl_start].to_vec();

    let template_r1 = r1_seq[tmpl_start..].to_vec();
    let template_r2 = r2_seq[tmpl_start..].to_vec();

    let umi_qual_r1 = r1_qual[..umi_len].to_vec();
    let umi_qual_r2 = r2_qual[..umi_len].to_vec();

    let qual_r1 = r1_qual[tmpl_start..].to_vec();
    let qual_r2 = r2_qual[tmpl_start..].to_vec();

    // ── 3. TemplateTooShort ───────────────────────────────────────────────
    if let Some(min_tmpl) = config.min_template_length {
        let t1 = template_r1.len();
        let t2 = template_r2.len();
        if t1 < min_tmpl || t2 < min_tmpl {
            let detail = format!("r1_template={t1}bp,r2_template={t2}bp,min={min_tmpl}");
            return ParseResult::Dropped {
                reason: DropReason::TemplateTooShort,
                detail,
            };
        }
    }

    // ── 4. LowUmiQuality ─────────────────────────────────────────────────
    // Quality bytes in FASTQ are Phred+33 encoded (ASCII).  Convert each byte
    // to a Phred score before comparing to the threshold.
    if let Some(min_q) = config.min_umi_quality {
        // Check R1 UMI qualities
        for (i, &raw) in umi_qual_r1.iter().enumerate() {
            let phred = raw.saturating_sub(33);
            if phred < min_q {
                let detail = format!("r1_umi_base={i},quality={phred},min={min_q}");
                return ParseResult::Dropped {
                    reason: DropReason::LowUmiQuality,
                    detail,
                };
            }
        }
        // Check R2 UMI qualities
        for (i, &raw) in umi_qual_r2.iter().enumerate() {
            let phred = raw.saturating_sub(33);
            if phred < min_q {
                let detail = format!("r2_umi_base={i},quality={phred},min={min_q}");
                return ParseResult::Dropped {
                    reason: DropReason::LowUmiQuality,
                    detail,
                };
            }
        }
    }

    // ── 5. Build canonical pair and strand ────────────────────────────────
    let canonical_umi = CanonicalUmiPair::new(umi_r1.clone(), umi_r2.clone());
    let strand = canonical_umi.strand_of_r1(&umi_r1);

    ParseResult::Ok(Box::new(ParsedReadPair {
        umi_r1,
        umi_r2,
        skip_r1,
        skip_r2,
        template_r1,
        template_r2,
        qual_r1,
        qual_r2,
        umi_qual_r1,
        umi_qual_r2,
        canonical_umi,
        strand,
    }))
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use kam_core::chemistry::ReadStructure;

    // Helpers -----------------------------------------------------------------

    fn default_config() -> ParserConfig {
        ParserConfig::default()
    }

    /// Build a read of the form [UMI][skip][template] with uniform quality `q`.
    fn make_read(umi: &[u8], skip: &[u8], template: &[u8]) -> (Vec<u8>, Vec<u8>) {
        let mut seq = Vec::new();
        seq.extend_from_slice(umi);
        seq.extend_from_slice(skip);
        seq.extend_from_slice(template);
        let qual = vec![b'I'; seq.len()]; // Phred 40
        (seq, qual)
    }

    // Test 1: Valid read pair with long template parses correctly — all fields
    // at right positions.
    #[test]
    fn valid_read_pair_long_template() {
        let (r1_seq, r1_qual) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNNNNNNNNNNNN");
        let (r2_seq, r2_qual) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNNNNNNNNNNNN");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config());
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.umi_r1, b"ACGTA");
                assert_eq!(p.umi_r2, b"TGCAT");
                assert_eq!(p.skip_r1.as_slice(), b"TG");
                assert_eq!(p.skip_r2.as_slice(), b"AC");
                assert_eq!(p.template_r1.len(), 20);
                assert_eq!(p.template_r2.len(), 20);
                assert_eq!(p.qual_r1.len(), 20);
                assert_eq!(p.qual_r2.len(), 20);
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Test 2: Read with 0 bp template (exactly umi+skip length) parses OK
    // when no min_template_length is set.
    #[test]
    fn zero_length_template_no_min_length() {
        let (r1_seq, r1_qual) = make_read(b"ACGTA", b"TG", b"");
        let (r2_seq, r2_qual) = make_read(b"TGCAT", b"AC", b"");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config());
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.template_r1.len(), 0);
                assert_eq!(p.template_r2.len(), 0);
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Test 3: Read shorter than umi+skip (e.g., 4 bp) returns
    // Dropped(ReadTooShort).
    #[test]
    fn read_shorter_than_umi_plus_skip_dropped() {
        let r1_seq = b"ACGT"; // 4 bp — too short (need 7)
        let r1_qual = b"IIII";
        let r2_seq = b"TGCATAGNNNNNNNNNN";
        let r2_qual = b"IIIIIIIIIIIIIIIII";
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &default_config());
        assert!(
            matches!(
                result,
                ParseResult::Dropped {
                    reason: DropReason::ReadTooShort,
                    ..
                }
            ),
            "Expected ReadTooShort"
        );
    }

    // Test 4: Read with 10 bp template returns Dropped(TemplateTooShort) when
    // min_template_length=20.
    #[test]
    fn template_too_short_when_min_set() {
        let (r1_seq, r1_qual) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNN"); // 10 bp tmpl
        let (r2_seq, r2_qual) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNN");
        let config = ParserConfig {
            min_template_length: Some(20),
            ..ParserConfig::default()
        };
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &config);
        assert!(
            matches!(
                result,
                ParseResult::Dropped {
                    reason: DropReason::TemplateTooShort,
                    ..
                }
            ),
            "Expected TemplateTooShort"
        );
    }

    // Test 5: Read with 10 bp template returns Ok when min_template_length is
    // None.
    #[test]
    fn template_10bp_ok_when_no_min_length() {
        let (r1_seq, r1_qual) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNN");
        let (r2_seq, r2_qual) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNN");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config());
        assert!(matches!(result, ParseResult::Ok(_)), "Expected Ok")
    }

    // Test 6: R1 and R2 of same length produce correct template lengths.
    #[test]
    fn r1_r2_same_length_correct_template_lengths() {
        let tmpl = b"ACGTACGTACGT"; // 12 bp
        let (r1_seq, r1_qual) = make_read(b"AAAAA", b"CC", tmpl);
        let (r2_seq, r2_qual) = make_read(b"TTTTT", b"GG", tmpl);
        assert_eq!(r1_seq.len(), r2_seq.len());
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config());
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.template_r1.len(), 12);
                assert_eq!(p.template_r2.len(), 12);
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Test 7: UMI extraction from positions 0..5.
    #[test]
    fn umi_extracted_from_positions_0_to_5() {
        let (r1_seq, r1_qual) = make_read(b"GGGGG", b"CC", b"AAAAAAAAAA");
        let (r2_seq, r2_qual) = make_read(b"CCCCC", b"TT", b"TTTTTTTTTT");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config());
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.umi_r1, b"GGGGG");
                assert_eq!(p.umi_r2, b"CCCCC");
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Test 8: Skip extraction from positions 5..7.
    #[test]
    fn skip_extracted_from_positions_5_to_7() {
        let (r1_seq, r1_qual) = make_read(b"AAAAA", b"XY", b"NNNNNNNNNN");
        let (r2_seq, r2_qual) = make_read(b"TTTTT", b"ZW", b"NNNNNNNNNN");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config());
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(&p.skip_r1, b"XY");
                assert_eq!(&p.skip_r2, b"ZW");
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Test 9: Template extraction starts at position 7.
    #[test]
    fn template_starts_at_position_7() {
        // Construct a read where we know exactly what the template should be.
        // Positions: 0..5 UMI="ACGTA", 5..7 skip="TG", 7..= template
        let r1_seq = b"ACGTATGCCCCCCCCCCC";
        let r1_qual = b"IIIIIIIIIIIIIIIIII";
        let r2_seq = b"TGCATACGGGGGGGGGGG";
        let r2_qual = b"IIIIIIIIIIIIIIIIII";
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &default_config());
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.template_r1, b"CCCCCCCCCCC");
                assert_eq!(p.template_r2, b"GGGGGGGGGGG");
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Test 10: Quality arrays match template length (not full read length).
    #[test]
    fn quality_arrays_match_template_length() {
        let tmpl = b"ACGTACGTACGT"; // 12 bp template
        let (r1_seq, r1_qual) = make_read(b"AAAAA", b"CC", tmpl);
        let (r2_seq, r2_qual) = make_read(b"TTTTT", b"GG", tmpl);
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config());
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.qual_r1.len(), p.template_r1.len());
                assert_eq!(p.qual_r2.len(), p.template_r2.len());
                // Neither should equal the full read length (19 bp each)
                assert_ne!(p.qual_r1.len(), r1_seq.len());
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Test 11: Canonical UMI and strand determined correctly.
    #[test]
    fn canonical_umi_and_strand_correct() {
        // "ACGTA" < "TGCAT" so umi_a = "ACGTA", R1 has the smaller UMI →
        // Forward strand.
        let (r1_seq, r1_qual) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNN");
        let (r2_seq, r2_qual) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNN");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config());
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.canonical_umi.umi_a, b"ACGTA");
                assert_eq!(p.canonical_umi.umi_b, b"TGCAT");
                assert_eq!(p.strand, Strand::Forward);
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Test 11b: Reverse strand when R1 UMI is the larger one.
    #[test]
    fn canonical_umi_reverse_strand_when_r1_umi_larger() {
        // "TGCAT" > "ACGTA" so canonical umi_a = "ACGTA", R1 has the larger →
        // Reverse strand.
        let (r1_seq, r1_qual) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNN");
        let (r2_seq, r2_qual) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNN");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config());
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.canonical_umi.umi_a, b"ACGTA");
                assert_eq!(p.strand, Strand::Reverse);
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Test 12: Low UMI quality base returns Dropped(LowUmiQuality) when
    // threshold set.
    #[test]
    fn low_umi_quality_dropped_when_threshold_set() {
        let r1_seq = b"ACGTATGNNNNNNNNNN";
        // Position 2 of the UMI has quality '#' = Phred 2
        let r1_qual = b"II#IIIIIIIIIIIIII";
        let r2_seq = b"TGCATAGNNNNNNNNNN";
        let r2_qual = b"IIIIIIIIIIIIIIIII";
        let config = ParserConfig {
            min_umi_quality: Some(30), // '#' = 2 < 30
            ..ParserConfig::default()
        };
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &config);
        assert!(
            matches!(
                result,
                ParseResult::Dropped {
                    reason: DropReason::LowUmiQuality,
                    ..
                }
            ),
            "Expected LowUmiQuality drop"
        );
    }

    // Test 13: Low UMI quality base returns Ok when no threshold set.
    #[test]
    fn low_umi_quality_ok_when_no_threshold() {
        let r1_seq = b"ACGTATGNNNNNNNNNN";
        // Position 2 of the UMI has quality '#' = Phred 2
        let r1_qual = b"II#IIIIIIIIIIIIII";
        let r2_seq = b"TGCATAGNNNNNNNNNN";
        let r2_qual = b"IIIIIIIIIIIIIIIII";
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &default_config());
        assert!(matches!(result, ParseResult::Ok(_)), "Expected Ok");
    }

    // Test 14: Detail string in Dropped contains useful info.
    #[test]
    fn dropped_detail_contains_useful_info() {
        // TemplateTooShort detail should mention template lengths and minimum.
        let (r1_seq, r1_qual) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNN"); // 10 bp
        let (r2_seq, r2_qual) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNN");
        let config = ParserConfig {
            min_template_length: Some(20),
            ..ParserConfig::default()
        };
        match parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &config) {
            ParseResult::Dropped {
                reason: DropReason::TemplateTooShort,
                detail,
            } => {
                assert!(
                    detail.contains("10"),
                    "detail should mention template length: {detail}"
                );
                assert!(
                    detail.contains("20"),
                    "detail should mention minimum: {detail}"
                );
            }
            other => panic!("Expected TemplateTooShort, got: {other:?}"),
        }

        // LowUmiQuality detail should mention which read/base and the quality.
        let r1_seq2 = b"ACGTATGNNNNNNNNNN";
        let r1_qual2 = b"II#IIIIIIIIIIIIII"; // '#' = Phred 2 at position 2
        let r2_seq2 = b"TGCATAGNNNNNNNNNN";
        let r2_qual2 = b"IIIIIIIIIIIIIIIII";
        let config2 = ParserConfig {
            min_umi_quality: Some(30),
            ..ParserConfig::default()
        };
        match parse_read_pair(r1_seq2, r1_qual2, r2_seq2, r2_qual2, &config2) {
            ParseResult::Dropped {
                reason: DropReason::LowUmiQuality,
                detail,
            } => {
                // detail should mention the quality value and the minimum
                assert!(
                    detail.contains("min=30"),
                    "detail should mention min quality: {detail}"
                );
            }
            other => panic!("Expected LowUmiQuality, got: {other:?}"),
        }

        // ReadTooShort detail should mention lengths.
        let short_seq = b"ACGT"; // 4 bp
        let short_qual = b"IIII";
        let long_seq = b"TGCATAGNNNNNNNNNN";
        let long_qual = b"IIIIIIIIIIIIIIIII";
        match parse_read_pair(
            short_seq,
            short_qual,
            long_seq,
            long_qual,
            &default_config(),
        ) {
            ParseResult::Dropped {
                reason: DropReason::ReadTooShort,
                detail,
            } => {
                assert!(
                    detail.contains('4'),
                    "detail should mention short read length: {detail}"
                );
                assert!(
                    detail.contains('7'),
                    "detail should mention min length: {detail}"
                );
            }
            other => panic!("Expected ReadTooShort, got: {other:?}"),
        }
    }

    // Test 15: Non-standard ReadStructure with umi_length != 5 parses
    // successfully — arbitrary UMI lengths are now supported.
    #[test]
    fn non_standard_umi_length_parses_ok() {
        let config = ParserConfig {
            read_structure: ReadStructure {
                umi_length: 8,
                skip_length: 2,
            },
            ..ParserConfig::default()
        };
        // 8 bp UMI + 2 bp skip + 11 bp template = 21 bp total
        let r1_seq = b"ACGTACGTATGNNNNNNNNNN";
        let r1_qual = b"IIIIIIIIIIIIIIIIIIIII";
        let r2_seq = b"TGCATGCAACGNNNNNNNNNN";
        let r2_qual = b"IIIIIIIIIIIIIIIIIIIII";
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &config);
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.umi_r1.len(), 8);
                assert_eq!(p.umi_r2.len(), 8);
                assert_eq!(p.umi_r1, b"ACGTACGT");
                assert_eq!(p.umi_r2, b"TGCATGCA");
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Test 16: Non-standard skip length (4 bp) parses correctly.
    // Verifies that skip lengths other than 2 bp are accepted.
    #[test]
    fn non_standard_skip_length_4bp_parses_ok() {
        let config = ParserConfig {
            read_structure: ReadStructure {
                umi_length: 5,
                skip_length: 4,
            },
            ..ParserConfig::default()
        };
        // 5 bp UMI + 4 bp skip + 9 bp template = 18 bp total.
        let r1_seq = b"ACGTATGCANNNNNNNN";
        let r1_qual = b"IIIIIIIIIIIIIIIII";
        let r2_seq = b"TGCATAGCTNNNNNNNNN";
        let r2_qual = b"IIIIIIIIIIIIIIIIII";
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &config);
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.umi_r1, b"ACGTA");
                assert_eq!(p.skip_r1.as_slice(), b"TGCA");
                assert_eq!(p.skip_r1.len(), 4);
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Test 16b: simplex_12bp chemistry (12 bp UMI, 0 bp skip) round-trips
    // through the parser without error.
    #[test]
    fn simplex_12bp_chemistry_roundtrip() {
        use kam_core::chemistry::ReadStructure;
        let config = ParserConfig {
            read_structure: ReadStructure::simplex_12bp(),
            ..ParserConfig::default()
        };
        // 12 bp UMI + 0 bp skip + 10 bp template = 22 bp total.
        let r1_seq = b"ACGTACGTACGTNNNNNNNNNN";
        let r1_qual = b"IIIIIIIIIIIIIIIIIIIIII";
        let r2_seq = b"TGCATGCATGCANNNNNNNNNN";
        let r2_qual = b"IIIIIIIIIIIIIIIIIIIIII";
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &config);
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.umi_r1, b"ACGTACGTACGT");
                assert_eq!(p.umi_r2, b"TGCATGCATGCA");
                assert_eq!(p.skip_r1.len(), 0, "simplex_12bp has 0 bp skip");
                assert_eq!(p.skip_r2.len(), 0, "simplex_12bp has 0 bp skip");
                assert_eq!(p.template_r1, b"NNNNNNNNNN");
                assert_eq!(p.template_r2, b"NNNNNNNNNN");
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Test 17: 12 bp UMI parses correctly end-to-end.
    #[test]
    fn twelve_bp_umi_parses_correctly() {
        let config = ParserConfig {
            read_structure: ReadStructure {
                umi_length: 12,
                skip_length: 2,
            },
            ..ParserConfig::default()
        };
        // 12 bp UMI + 2 bp skip + 10 bp template = 24 bp total
        let r1_seq = b"ACGTACGTACGTTGNNNNNNNNNN";
        let r1_qual = b"IIIIIIIIIIIIIIIIIIIIIIII";
        let r2_seq = b"TGCATGCATGCAACNNNNNNNNNN";
        let r2_qual = b"IIIIIIIIIIIIIIIIIIIIIIII";
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &config);
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.umi_r1.len(), 12);
                assert_eq!(p.umi_r2.len(), 12);
                assert_eq!(p.umi_r1, b"ACGTACGTACGT");
                assert_eq!(p.umi_r2, b"TGCATGCATGCA");
                // "ACGTACGTACGT" < "TGCATGCATGCA" — forward strand
                assert_eq!(p.strand, Strand::Forward);
                assert_eq!(p.canonical_umi.umi_a, b"ACGTACGTACGT");
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Additional: ParseStats default values are all zero.
    #[test]
    fn parse_stats_default_all_zero() {
        let stats = ParseStats::default();
        assert_eq!(stats.n_processed, 0);
        assert_eq!(stats.n_passed, 0);
        assert_eq!(stats.n_read_too_short, 0);
        assert_eq!(stats.n_template_too_short, 0);
        assert_eq!(stats.n_low_umi_quality, 0);
    }

    // Additional: default ParserConfig uses the Twist UMI duplex preset.
    #[test]
    fn default_config_has_twist_preset() {
        let config = ParserConfig::default();
        assert_eq!(config.read_structure.umi_length, 5);
        assert_eq!(config.read_structure.skip_length, 2);
    }

    // ── Edge-case tests ─────────────────────────────────────────────────────

    // UMI base exactly at quality threshold passes; one below fails.
    // Phred+33 encoding: Q=20 is ASCII 53 = '5'.
    #[test]
    fn umi_quality_at_threshold_passes() {
        let config = ParserConfig {
            min_umi_quality: Some(20),
            ..ParserConfig::default()
        };
        let q20_char = (20 + 33) as u8; // ASCII 53 = '5'
        let (r1_seq, _) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNN");
        let mut r1_qual = vec![b'I'; r1_seq.len()];
        for q in r1_qual.iter_mut().take(5) {
            *q = q20_char;
        }
        let (r2_seq, r2_qual) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNN");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &config);
        assert!(
            matches!(result, ParseResult::Ok(_)),
            "UMI base at exactly threshold should pass"
        );
    }

    #[test]
    fn umi_quality_one_below_threshold_fails() {
        let config = ParserConfig {
            min_umi_quality: Some(20),
            ..ParserConfig::default()
        };
        let q19_char = (19 + 33) as u8;
        let (r1_seq, _) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNN");
        let mut r1_qual = vec![b'I'; r1_seq.len()];
        r1_qual[0] = q19_char; // First UMI base below threshold.
        let (r2_seq, r2_qual) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNN");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &config);
        assert!(
            matches!(
                result,
                ParseResult::Dropped {
                    reason: DropReason::LowUmiQuality,
                    ..
                }
            ),
            "UMI base one below threshold should be dropped"
        );
    }

    // All UMI bases fail quality: molecule is dropped with correct reason.
    #[test]
    fn all_umi_bases_fail_quality_dropped() {
        let config = ParserConfig {
            min_umi_quality: Some(30),
            ..ParserConfig::default()
        };
        let (r1_seq, _) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNN");
        let mut r1_qual = vec![b'I'; r1_seq.len()];
        for q in r1_qual.iter_mut().take(5) {
            *q = b'#'; // Q=2, well below threshold.
        }
        let (r2_seq, r2_qual) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNN");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &config);
        match result {
            ParseResult::Dropped {
                reason: DropReason::LowUmiQuality,
                detail,
            } => {
                assert!(
                    detail.contains("r1_umi_base=0"),
                    "detail should report base index 0: {detail}"
                );
            }
            other => panic!("Expected LowUmiQuality drop, got: {other:?}"),
        }
    }

    // Skip length 1 bp: parses correctly (variable skip is now supported).
    #[test]
    fn skip_length_1bp_parses_ok() {
        let config = ParserConfig {
            read_structure: ReadStructure {
                umi_length: 5,
                skip_length: 1,
            },
            ..ParserConfig::default()
        };
        // 5 bp UMI + 1 bp skip + 10 bp template = 16 bp total.
        let r1_seq = b"ACGTAXNNNNNNNNNN";
        let r1_qual = b"IIIIIIIIIIIIIIII";
        let r2_seq = b"TGCATANNNNNNNNNN";
        let r2_qual = b"IIIIIIIIIIIIIIII";
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &config);
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.skip_r1.as_slice(), b"X");
                assert_eq!(p.skip_r2.as_slice(), b"A");
                assert_eq!(p.skip_r1.len(), 1);
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Skip length 3 bp: parses correctly (variable skip is now supported).
    #[test]
    fn skip_length_3bp_parses_ok() {
        let config = ParserConfig {
            read_structure: ReadStructure {
                umi_length: 5,
                skip_length: 3,
            },
            ..ParserConfig::default()
        };
        // 5 bp UMI + 3 bp skip + 10 bp template = 18 bp total.
        let r1_seq = b"ACGTATGCNNNNNNNNNN";
        let r1_qual = b"IIIIIIIIIIIIIIIIII";
        let r2_seq = b"TGCATACGNNNNNNNNNNN";
        let r2_qual = b"IIIIIIIIIIIIIIIIIII";
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &config);
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.skip_r1.as_slice(), b"TGC");
                assert_eq!(p.skip_r1.len(), 3);
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Template exactly at minimum length: passes.
    #[test]
    fn template_exactly_at_min_length_passes() {
        let config = ParserConfig {
            min_template_length: Some(10),
            ..ParserConfig::default()
        };
        let (r1_seq, r1_qual) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNN"); // 10 bp template
        let (r2_seq, r2_qual) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNN");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &config);
        assert!(
            matches!(result, ParseResult::Ok(_)),
            "template at exactly min_template_length should pass"
        );
    }

    // Template one base shorter than minimum: dropped.
    #[test]
    fn template_one_below_min_length_dropped() {
        let config = ParserConfig {
            min_template_length: Some(10),
            ..ParserConfig::default()
        };
        let (r1_seq, r1_qual) = make_read(b"ACGTA", b"TG", b"NNNNNNNNN"); // 9 bp template
        let (r2_seq, r2_qual) = make_read(b"TGCAT", b"AC", b"NNNNNNNNN");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &config);
        assert!(
            matches!(
                result,
                ParseResult::Dropped {
                    reason: DropReason::TemplateTooShort,
                    ..
                }
            ),
            "template one below minimum should be dropped"
        );
    }

    // Canonical UMI pair is symmetric: parse(A+B) and parse(B+A) give same canonical.
    #[test]
    fn canonical_umi_pair_is_symmetric() {
        let (r1_a, q1_a) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNN");
        let (r2_a, q2_a) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNN");
        let result_a = parse_read_pair(&r1_a, &q1_a, &r2_a, &q2_a, &default_config());

        let (r1_b, q1_b) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNN");
        let (r2_b, q2_b) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNN");
        let result_b = parse_read_pair(&r1_b, &q1_b, &r2_b, &q2_b, &default_config());

        let canon_a = match result_a {
            ParseResult::Ok(p) => p.canonical_umi.clone(),
            ParseResult::Dropped { .. } => panic!("expected Ok"),
        };
        let canon_b = match result_b {
            ParseResult::Ok(p) => p.canonical_umi.clone(),
            ParseResult::Dropped { .. } => panic!("expected Ok"),
        };

        assert_eq!(
            canon_a, canon_b,
            "canonical UMI pair must be the same regardless of R1/R2 order"
        );
    }

    // UMI with N bases: parses successfully when no quality threshold is set.
    #[test]
    fn umi_with_n_bases_parses_when_no_quality_threshold() {
        let (r1_seq, r1_qual) = make_read(b"NNNNN", b"TG", b"NNNNNNNNNN");
        let (r2_seq, r2_qual) = make_read(b"NNNNN", b"AC", b"NNNNNNNNNN");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config());
        assert!(
            matches!(result, ParseResult::Ok(_)),
            "N bases in UMI should parse successfully without quality filter"
        );
    }

    // All-N template: parses successfully.
    #[test]
    fn all_n_template_parses_successfully() {
        let (r1_seq, r1_qual) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNN");
        let (r2_seq, r2_qual) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNN");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config());
        match result {
            ParseResult::Ok(p) => {
                assert!(
                    p.template_r1.iter().all(|&b| b == b'N'),
                    "template should be all N"
                );
            }
            ParseResult::Dropped { .. } => panic!("expected Ok for all-N template"),
        }
    }

    // Read shorter than UMI + skip: dropped, not a panic.
    #[test]
    fn read_exactly_one_byte_short_dropped_not_panic() {
        let r1_seq = b"ACGTAC"; // 6 bp, one short of 7
        let r1_qual = b"IIIIII";
        let r2_seq = b"TGCATAGNNNNNNNNNN";
        let r2_qual = b"IIIIIIIIIIIIIIIII";
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &default_config());
        assert!(
            matches!(
                result,
                ParseResult::Dropped {
                    reason: DropReason::ReadTooShort,
                    ..
                }
            ),
            "6 bp read (one short of 7) should be dropped as ReadTooShort"
        );
    }

    // R1 and R2 with different template lengths: both parse successfully.
    #[test]
    fn r1_r2_different_template_lengths_both_parse() {
        let (r1_seq, r1_qual) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNN"); // 10 bp template
        let (r2_seq, r2_qual) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNNNNNNN"); // 15 bp template
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config());
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.template_r1.len(), 10, "R1 template should be 10 bp");
                assert_eq!(p.template_r2.len(), 15, "R2 template should be 15 bp");
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // R2 UMI quality failure is detected even when R1 UMI quality is fine.
    #[test]
    fn r2_umi_quality_failure_detected() {
        let config = ParserConfig {
            min_umi_quality: Some(30),
            ..ParserConfig::default()
        };
        let (r1_seq, r1_qual) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNN");
        let (r2_seq, _) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNN");
        let mut r2_qual = vec![b'I'; r2_seq.len()];
        r2_qual[2] = b'#'; // Third UMI base at Q=2, below threshold.
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &config);
        match result {
            ParseResult::Dropped {
                reason: DropReason::LowUmiQuality,
                detail,
            } => {
                assert!(
                    detail.contains("r2_umi_base"),
                    "detail should mention R2: {detail}"
                );
            }
            other => panic!("Expected LowUmiQuality for R2, got: {other:?}"),
        }
    }

    // Read exactly at umi+skip length (0 bp template) with min_template_length
    // set: dropped as TemplateTooShort, not ReadTooShort.
    #[test]
    fn zero_template_with_min_template_length_dropped_as_template_too_short() {
        let config = ParserConfig {
            min_template_length: Some(1),
            ..ParserConfig::default()
        };
        let (r1_seq, r1_qual) = make_read(b"ACGTA", b"TG", b""); // 0 bp template
        let (r2_seq, r2_qual) = make_read(b"TGCAT", b"AC", b"");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &config);
        assert!(
            matches!(
                result,
                ParseResult::Dropped {
                    reason: DropReason::TemplateTooShort,
                    ..
                }
            ),
            "0 bp template with min=1 should be TemplateTooShort, not ReadTooShort"
        );
    }

    // Empty read (0 bytes) is dropped as ReadTooShort, not a panic.
    #[test]
    fn empty_read_dropped_not_panic() {
        let result = parse_read_pair(
            b"",
            b"",
            b"TGCATAGNNNNNNNNNN",
            b"IIIIIIIIIIIIIIIII",
            &default_config(),
        );
        assert!(
            matches!(
                result,
                ParseResult::Dropped {
                    reason: DropReason::ReadTooShort,
                    ..
                }
            ),
            "empty read should be dropped as ReadTooShort"
        );
    }

    // Both R1 and R2 empty: dropped as ReadTooShort.
    #[test]
    fn both_reads_empty_dropped() {
        let result = parse_read_pair(b"", b"", b"", b"", &default_config());
        assert!(
            matches!(
                result,
                ParseResult::Dropped {
                    reason: DropReason::ReadTooShort,
                    ..
                }
            ),
            "both reads empty should be dropped as ReadTooShort"
        );
    }

    // UMI quality check order: R1 is checked before R2. When both have low
    // quality, the detail mentions R1.
    #[test]
    fn umi_quality_check_order_r1_before_r2() {
        let config = ParserConfig {
            min_umi_quality: Some(30),
            ..ParserConfig::default()
        };
        let (r1_seq, _) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNN");
        let mut r1_qual = vec![b'#'; r1_seq.len()]; // All Q=2 for R1
        for q in r1_qual.iter_mut().skip(7) {
            *q = b'I';
        }
        let (r2_seq, _) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNN");
        let mut r2_qual = vec![b'#'; r2_seq.len()];
        for q in r2_qual.iter_mut().skip(7) {
            *q = b'I';
        }
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &config);
        match result {
            ParseResult::Dropped {
                reason: DropReason::LowUmiQuality,
                detail,
            } => {
                assert!(
                    detail.contains("r1_umi_base"),
                    "should detect R1 failure first: {detail}"
                );
            }
            other => panic!("Expected LowUmiQuality, got: {other:?}"),
        }
    }

    // 9 bp UMI with non-standard ReadStructure parses correctly.
    #[test]
    fn nine_bp_umi_parses_correctly() {
        let config = ParserConfig {
            read_structure: ReadStructure {
                umi_length: 9,
                skip_length: 2,
            },
            ..ParserConfig::default()
        };
        // 9 bp UMI + 2 bp skip + 10 bp template = 21 bp total
        let r1_seq = b"ACGTACGTAATNNNNNNNNNN";
        let r1_qual = b"IIIIIIIIIIIIIIIIIIIII";
        let r2_seq = b"TGCATGCATATNNNNNNNNNN";
        let r2_qual = b"IIIIIIIIIIIIIIIIIIIII";
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &config);
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.umi_r1.len(), 9);
                assert_eq!(p.umi_r2.len(), 9);
                assert_eq!(p.umi_r1, b"ACGTACGTA");
                assert_eq!(p.umi_r2, b"TGCATGCAT");
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // UMI quality fields are correctly extracted and have the right length.
    #[test]
    fn umi_quality_fields_have_correct_length_and_values() {
        let (r1_seq, _) = make_read(b"ACGTA", b"TG", b"NNNNNNNNNN");
        let mut r1_qual = vec![b'I'; r1_seq.len()];
        r1_qual[0] = b'5'; // Q=20
        r1_qual[1] = b'?'; // Q=30
        r1_qual[2] = b'I'; // Q=40
        r1_qual[3] = b'5'; // Q=20
        r1_qual[4] = b'?'; // Q=30
        let (r2_seq, r2_qual) = make_read(b"TGCAT", b"AC", b"NNNNNNNNNN");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config());
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.umi_qual_r1.len(), 5, "UMI quality should be 5 bytes");
                assert_eq!(p.umi_qual_r1[0], b'5');
                assert_eq!(p.umi_qual_r1[1], b'?');
                assert_eq!(p.umi_qual_r2.len(), 5);
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Identical UMIs on R1 and R2: canonical pair has umi_a == umi_b, strand
    // is Forward (since R1 UMI == umi_a when they are equal).
    #[test]
    fn identical_umi_r1_r2_canonical_pair() {
        let (r1_seq, r1_qual) = make_read(b"AAAAA", b"TG", b"NNNNNNNNNN");
        let (r2_seq, r2_qual) = make_read(b"AAAAA", b"AC", b"NNNNNNNNNN");
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config());
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.canonical_umi.umi_a, b"AAAAA");
                assert_eq!(p.canonical_umi.umi_b, b"AAAAA");
                assert_eq!(
                    p.strand,
                    Strand::Forward,
                    "R1 UMI == umi_a when equal, so strand should be Forward"
                );
            }
            ParseResult::Dropped { reason, detail } => {
                panic!("Expected Ok, got Dropped({reason:?}): {detail}");
            }
        }
    }

    // Filter priority: ReadTooShort is checked before TemplateTooShort and
    // LowUmiQuality.
    #[test]
    fn filter_priority_read_too_short_first() {
        let config = ParserConfig {
            min_template_length: Some(100),
            min_umi_quality: Some(40),
            ..ParserConfig::default()
        };
        let result = parse_read_pair(b"ACGT", b"IIII", b"ACGT", b"IIII", &config);
        assert!(
            matches!(
                result,
                ParseResult::Dropped {
                    reason: DropReason::ReadTooShort,
                    ..
                }
            ),
            "ReadTooShort should take priority over other filters"
        );
    }
}
