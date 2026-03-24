//! Twist UMI FASTQ read pair parser.
//!
//! Extracts UMI, skip bases, and template from raw R1/R2 reads following the
//! Twist UMI duplex `5M2S+T` read structure.  The parser is k-agnostic: it
//! does **not** enforce a minimum template length by default.  Short templates
//! are parsed successfully; downstream stages (e.g. k-mer extraction) decide
//! whether to use them.

use kam_core::chemistry::ReadStructure;
use kam_core::molecule::{CanonicalUmiPair, Strand};

// ── Error type ────────────────────────────────────────────────────────────────

/// Errors returned by [`parse_read_pair`].
///
/// Currently this enum is uninhabited: the parser stores UMI and skip bytes as
/// `Vec<u8>`, so any `ReadStructure` length is accepted without error.  The
/// type is preserved so that callers need not change their `Result` handling if
/// future validation is added.
///
/// # Example
/// ```
/// use kam_assemble::parser::{ParserConfig, ParseResult, parse_read_pair};
/// use kam_core::chemistry::ReadStructure;
///
/// // Non-standard read structure with UMI length 8 — accepted without error.
/// let config = ParserConfig {
///     read_structure: ReadStructure { umi_length: 8, skip_length: 2 },
///     ..ParserConfig::default()
/// };
/// let r1_seq  = b"ACGTACGTATGNNNNNNNNNN";
/// let r1_qual = b"IIIIIIIIIIIIIIIIIIIII";
/// let r2_seq  = b"TGCATACGTAGNNNNNNNNNN";
/// let r2_qual = b"IIIIIIIIIIIIIIIIIIIII";
/// let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &config);
/// assert!(matches!(result, Ok(ParseResult::Ok(_))));
/// ```
#[derive(Debug, thiserror::Error)]
pub enum ParseError {}

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
/// `ParsedReadPair` is boxed to keep the enum size small; `Dropped` is the
/// hot-path rejection case and should remain unboxed.
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
/// UMI and skip bytes are stored as `Vec<u8>`, so any `ReadStructure` length
/// is accepted.
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
/// # Errors
///
/// Returns `Err` only when a future validation is added to [`ParseError`].
/// Currently this function is infallible (returns `Ok` for all valid inputs).
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
/// let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &config).unwrap();
/// assert!(matches!(result, ParseResult::Ok(_)));
/// ```
pub fn parse_read_pair(
    r1_seq: &[u8],
    r1_qual: &[u8],
    r2_seq: &[u8],
    r2_qual: &[u8],
    config: &ParserConfig,
) -> Result<ParseResult, ParseError> {
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
        return Ok(ParseResult::Dropped {
            reason: DropReason::ReadTooShort,
            detail,
        });
    }

    // ── 2. Extract UMI, skip, template ───────────────────────────────────
    // The ReadTooShort guard above ensures both reads are at least `tmpl_start`
    // bytes long, so all slice bounds are valid.  UMI and skip are stored as
    // Vec<u8> to support arbitrary lengths.
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
            return Ok(ParseResult::Dropped {
                reason: DropReason::TemplateTooShort,
                detail,
            });
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
                return Ok(ParseResult::Dropped {
                    reason: DropReason::LowUmiQuality,
                    detail,
                });
            }
        }
        // Check R2 UMI qualities
        for (i, &raw) in umi_qual_r2.iter().enumerate() {
            let phred = raw.saturating_sub(33);
            if phred < min_q {
                let detail = format!("r2_umi_base={i},quality={phred},min={min_q}");
                return Ok(ParseResult::Dropped {
                    reason: DropReason::LowUmiQuality,
                    detail,
                });
            }
        }
    }

    // ── 5. Build canonical pair and strand ────────────────────────────────
    let canonical_umi = CanonicalUmiPair::new(umi_r1.clone(), umi_r2.clone());
    let strand = canonical_umi.strand_of_r1(&umi_r1);

    Ok(ParseResult::Ok(Box::new(ParsedReadPair {
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
    })))
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

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
        let result =
            parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config()).unwrap();
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.umi_r1, b"ACGTA");
                assert_eq!(p.umi_r2, b"TGCAT");
                assert_eq!(p.skip_r1, b"TG");
                assert_eq!(p.skip_r2, b"AC");
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
        let result =
            parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config()).unwrap();
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
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &default_config()).unwrap();
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
        let result = parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &config).unwrap();
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
        let result =
            parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config()).unwrap();
        assert!(matches!(result, ParseResult::Ok(_)), "Expected Ok");
    }

    // Test 6: R1 and R2 of same length produce correct template lengths.
    #[test]
    fn r1_r2_same_length_correct_template_lengths() {
        let tmpl = b"ACGTACGTACGT"; // 12 bp
        let (r1_seq, r1_qual) = make_read(b"AAAAA", b"CC", tmpl);
        let (r2_seq, r2_qual) = make_read(b"TTTTT", b"GG", tmpl);
        assert_eq!(r1_seq.len(), r2_seq.len());
        let result =
            parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config()).unwrap();
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
        let result =
            parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config()).unwrap();
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
        let result =
            parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config()).unwrap();
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.skip_r1, b"XY");
                assert_eq!(p.skip_r2, b"ZW");
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
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &default_config()).unwrap();
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
        let result =
            parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config()).unwrap();
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
        let result =
            parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config()).unwrap();
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
        let result =
            parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &default_config()).unwrap();
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
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &config).unwrap();
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
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &default_config()).unwrap();
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
        match parse_read_pair(&r1_seq, &r1_qual, &r2_seq, &r2_qual, &config).unwrap() {
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
        match parse_read_pair(r1_seq2, r1_qual2, r2_seq2, r2_qual2, &config2).unwrap() {
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
        )
        .unwrap()
        {
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

    // Test 15: Non-standard ReadStructure with umi_length != 5 is now accepted
    // — UMI is stored as Vec<u8>, so any length parses successfully.
    #[test]
    fn non_standard_umi_length_parses_ok() {
        use kam_core::chemistry::ReadStructure;
        let config = ParserConfig {
            read_structure: ReadStructure {
                umi_length: 8,
                skip_length: 2,
            },
            ..ParserConfig::default()
        };
        // 8 bp UMI + 2 bp skip + 11 bp template = 21 bp total.
        let r1_seq = b"ACGTACGTATGNNNNNNNNNN";
        let r1_qual = b"IIIIIIIIIIIIIIIIIIIII";
        let r2_seq = b"TGCATACGTAGNNNNNNNNNN";
        let r2_qual = b"IIIIIIIIIIIIIIIIIIIII";
        let result = parse_read_pair(r1_seq, r1_qual, r2_seq, r2_qual, &config).unwrap();
        match result {
            ParseResult::Ok(p) => {
                assert_eq!(p.umi_r1.len(), 8, "UMI should be 8 bp");
                assert_eq!(p.skip_r1.len(), 2, "skip should be 2 bp");
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
}
