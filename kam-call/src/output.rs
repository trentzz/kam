//! Writers for variant calls in TSV, CSV, JSON, and VCF formats.
//!
//! All writers accept a `&[VariantCall]` slice and a `&mut dyn Write` target,
//! enabling use with files, sockets, or in-memory buffers equally.

use std::io::{self, Write};

use serde_json::{json, Value};

use crate::allele::extract_minimal_allele;
use crate::caller::{VariantCall, VariantFilter, VariantType};
use crate::fusion::FusionCall;

// ─── Public types ─────────────────────────────────────────────────────────────

/// Supported output formats for variant calls.
///
/// # Example
/// ```
/// use kam_call::output::OutputFormat;
/// let fmt = OutputFormat::Tsv;
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    /// Tab-separated values (km-compatible).
    Tsv,
    /// Comma-separated values.
    Csv,
    /// JSON array of variant objects.
    Json,
    /// VCF 4.3 with custom INFO fields.
    Vcf,
}

// ─── Public functions ─────────────────────────────────────────────────────────

/// Write variant calls in any supported format.
///
/// Dispatches to the appropriate format-specific writer.
///
/// `ml_filter_threshold` is the probability threshold at or above which a call
/// is labelled `ML_PASS` in the `ml_filter` column. Pass `0.5` when no model
/// is loaded. Use `scorer.meta.ml_pass_threshold` when a model is loaded.
///
/// # Example
/// ```
/// use kam_call::output::{write_variants, OutputFormat};
/// use kam_call::caller::{VariantCall, VariantType, VariantFilter};
///
/// let calls: Vec<VariantCall> = vec![];
/// let mut buf = Vec::new();
/// write_variants(&calls, OutputFormat::Tsv, &mut buf, 0.5).unwrap();
/// assert!(!buf.is_empty());
/// ```
pub fn write_variants(
    calls: &[VariantCall],
    format: OutputFormat,
    writer: &mut dyn Write,
    ml_filter_threshold: f64,
) -> io::Result<()> {
    match format {
        OutputFormat::Tsv => write_delimited(calls, '\t', writer, ml_filter_threshold),
        OutputFormat::Csv => write_delimited(calls, ',', writer, ml_filter_threshold),
        OutputFormat::Json => write_json_with_threshold(calls, writer, ml_filter_threshold),
        OutputFormat::Vcf => write_vcf(calls, writer),
    }
}

/// Write variant calls as tab-separated values (km-compatible format).
///
/// Produces a header line followed by one data line per variant. Sequences
/// are encoded as UTF-8 strings; invalid bytes are replaced with `?`.
///
/// # Example
/// ```
/// use kam_call::output::write_tsv;
/// let calls = vec![];
/// let mut buf = Vec::new();
/// write_tsv(&calls, &mut buf).unwrap();
/// let text = String::from_utf8(buf).unwrap();
/// assert!(text.starts_with("target_id\t"));
/// ```
pub fn write_tsv(calls: &[VariantCall], writer: &mut dyn Write) -> io::Result<()> {
    write_delimited(calls, '\t', writer, 0.5)
}

/// Write variant calls as comma-separated values.
///
/// # Example
/// ```
/// use kam_call::output::write_csv;
/// let calls = vec![];
/// let mut buf = Vec::new();
/// write_csv(&calls, &mut buf).unwrap();
/// let text = String::from_utf8(buf).unwrap();
/// assert!(text.starts_with("target_id,"));
/// ```
pub fn write_csv(calls: &[VariantCall], writer: &mut dyn Write) -> io::Result<()> {
    write_delimited(calls, ',', writer, 0.5)
}

/// Write variant calls as a JSON array of objects.
///
/// Sequences are represented as ASCII strings. An empty slice produces `[]`.
///
/// # Example
/// ```
/// use kam_call::output::write_json;
/// let calls = vec![];
/// let mut buf = Vec::new();
/// write_json(&calls, &mut buf).unwrap();
/// let text = String::from_utf8(buf).unwrap();
/// assert_eq!(text.trim(), "[]");
/// ```
pub fn write_json(calls: &[VariantCall], writer: &mut dyn Write) -> io::Result<()> {
    write_json_with_threshold(calls, writer, 0.5)
}

/// Write variant calls as minimal VCF 4.3.
///
/// Uses `target_id` as CHROM and position 1 as POS (alignment-free).
/// Custom INFO fields carry all statistical evidence.
///
/// # Example
/// ```
/// use kam_call::output::write_vcf;
/// let calls = vec![];
/// let mut buf = Vec::new();
/// write_vcf(&calls, &mut buf).unwrap();
/// let text = String::from_utf8(buf).unwrap();
/// assert!(text.contains("##fileformat=VCF"));
/// ```
pub fn write_vcf(calls: &[VariantCall], writer: &mut dyn Write) -> io::Result<()> {
    // VCF header
    writeln!(writer, "##fileformat=VCFv4.3")?;
    writeln!(writer, "##source=kam-call")?;
    writeln!(
        writer,
        "##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Variant allele frequency\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=VAF_LO,Number=1,Type=Float,Description=\"VAF 95% CI lower\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=VAF_HI,Number=1,Type=Float,Description=\"VAF 95% CI upper\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=NREF,Number=1,Type=Integer,Description=\"Molecules supporting REF\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=NALT,Number=1,Type=Integer,Description=\"Molecules supporting ALT\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=NDUPALT,Number=1,Type=Integer,Description=\"Duplex molecules supporting ALT\">"
    )?;
    writeln!(writer, "##INFO=<ID=NSIMALT,Number=1,Type=Integer,Description=\"Simplex molecules supporting ALT\">")?;
    writeln!(
        writer,
        "##INFO=<ID=SBP,Number=1,Type=Float,Description=\"Strand bias Fisher p-value\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=CONF,Number=1,Type=Float,Description=\"Posterior confidence\">"
    )?;
    writeln!(
        writer,
        "##FILTER=<ID=PASS,Description=\"All filters passed\">"
    )?;
    writeln!(
        writer,
        "##FILTER=<ID=StrandBias,Description=\"Strand bias detected\">"
    )?;
    writeln!(
        writer,
        "##FILTER=<ID=LowConfidence,Description=\"Low posterior confidence\">"
    )?;
    writeln!(
        writer,
        "##FILTER=<ID=LowDuplex,Description=\"Insufficient duplex support\">"
    )?;
    writeln!(
        writer,
        "##FILTER=<ID=CollisionRisk,Description=\"UMI collision risk\">"
    )?;
    writeln!(
        writer,
        "##FILTER=<ID=HighVaf,Description=\"VAF exceeds maximum threshold (likely germline)\">"
    )?;
    writeln!(
        writer,
        "##FILTER=<ID=NotTargeted,Description=\"Allele not in --target-variants set (tumour-informed monitoring mode)\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"
    )?;
    writeln!(writer, "##ALT=<ID=DEL,Description=\"Deletion\">")?;
    writeln!(writer, "##ALT=<ID=DUP,Description=\"Tandem duplication\">")?;
    writeln!(writer, "##ALT=<ID=INV,Description=\"Inversion\">")?;
    writeln!(
        writer,
        "##ALT=<ID=BND,Description=\"Breakend (fusion junction)\">"
    )?;
    writeln!(
        writer,
        "##ALT=<ID=INVDEL,Description=\"Inversion with flanking deletion\">"
    )?;
    writeln!(
        writer,
        "##ALT=<ID=INS,Description=\"Novel insertion (not a tandem duplication)\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=MATEID,Number=1,Type=String,Description=\"ID of the mate breakend\">"
    )?;
    writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;

    for call in calls {
        let filter_str = filter_to_str(call.filter);
        let base_info = format!(
            "VAF={:.6};VAF_LO={:.6};VAF_HI={:.6};NREF={};NALT={};NDUPALT={};NSIMALT={};SBP={:.6};CONF={:.6}",
            call.vaf,
            call.vaf_ci_low,
            call.vaf_ci_high,
            call.n_molecules_ref,
            call.n_molecules_alt,
            call.n_duplex_alt,
            call.n_simplex_alt,
            call.strand_bias_p,
            call.confidence,
        );

        // Structural variants use symbolic allele notation in VCF.
        if is_sv_type(call.variant_type) {
            write_sv_vcf_record(call, filter_str, &base_info, writer)?;
            continue;
        }

        // When target_id encodes a genomic coordinate (chrN:START-END), emit a
        // properly placed VCF record with minimal left-normalised alleles.
        // Fall back to the alignment-free format (target_id as CHROM, POS=1)
        // when the target_id cannot be parsed.
        if let Some(a) =
            extract_minimal_allele(&call.target_id, &call.ref_sequence, &call.alt_sequence)
        {
            writeln!(
                writer,
                "{}\t{}\t.\t{}\t{}\t.\t{}\t{}",
                a.chrom, a.pos, a.ref_allele, a.alt_allele, filter_str, base_info
            )?;
        } else {
            let ref_str = bytes_to_str(&call.ref_sequence);
            let alt_str = bytes_to_str(&call.alt_sequence);
            writeln!(
                writer,
                "{}\t1\t.\t{}\t{}\t.\t{}\t{}",
                call.target_id, ref_str, alt_str, filter_str, base_info
            )?;
        }
    }
    Ok(())
}

/// Return true when the variant type is a structural variant.
fn is_sv_type(vt: VariantType) -> bool {
    matches!(
        vt,
        VariantType::LargeDeletion
            | VariantType::TandemDuplication
            | VariantType::Inversion
            | VariantType::Fusion
            | VariantType::InvDel
            | VariantType::NovelInsertion
    )
}

/// Write a VCF record for a structural variant using symbolic allele notation.
///
/// Uses the minimal allele extractor to get the correct anchor coordinate, then
/// emits `<DEL>`, `<DUP>`, or `<INV>` in the ALT field with SVTYPE and SVLEN
/// appended to the INFO string.
fn write_sv_vcf_record(
    call: &VariantCall,
    filter_str: &str,
    base_info: &str,
    writer: &mut dyn Write,
) -> io::Result<()> {
    let (sym_alt, svtype_str) = match call.variant_type {
        VariantType::LargeDeletion => ("<DEL>", "DEL"),
        VariantType::TandemDuplication => ("<DUP>", "DUP"),
        VariantType::Inversion => ("<INV>", "INV"),
        VariantType::Fusion => ("<BND>", "BND"),
        VariantType::InvDel => ("<INVDEL>", "INVDEL"),
        VariantType::NovelInsertion => ("<INS>", "INS"),
        _ => unreachable!("write_sv_vcf_record called with non-SV type"),
    };

    // SVLEN: negative for deletions, positive for insertions, 0 for inversions and fusions.
    let svlen: i64 = match call.variant_type {
        VariantType::LargeDeletion => {
            -((call.ref_sequence.len() as i64) - (call.alt_sequence.len() as i64))
        }
        VariantType::TandemDuplication | VariantType::NovelInsertion => {
            (call.alt_sequence.len() as i64) - (call.ref_sequence.len() as i64)
        }
        VariantType::InvDel => {
            -((call.ref_sequence.len() as i64) - (call.alt_sequence.len() as i64))
        }
        VariantType::Inversion | VariantType::Fusion => 0,
        _ => unreachable!(),
    };

    let info = format!("{base_info};SVTYPE={svtype_str};SVLEN={svlen}");

    // Use the minimal allele extractor to get the correct coordinate; fall back to POS=1.
    if let Some(a) = extract_minimal_allele(&call.target_id, &call.ref_sequence, &call.alt_sequence)
    {
        // For symbolic alleles the REF field must contain only the anchor base.
        let ref_anchor = &a.ref_allele[..1];
        writeln!(
            writer,
            "{}\t{}\t.\t{}\t{}\t.\t{}\t{}",
            a.chrom, a.pos, ref_anchor, sym_alt, filter_str, info
        )?;
    } else {
        // Alignment-free fallback: REF = first base of ref_sequence.
        let ref_anchor = if call.ref_sequence.is_empty() {
            "N".to_string()
        } else {
            (call.ref_sequence[0] as char).to_string()
        };
        writeln!(
            writer,
            "{}\t1\t.\t{}\t{}\t.\t{}\t{}",
            call.target_id, ref_anchor, sym_alt, filter_str, info
        )?;
    }
    Ok(())
}

/// Write two paired VCF BND records for a fusion call.
///
/// Emits one record per partner breakpoint following VCF 4.3 BND notation.
/// The two records are linked by `MATEID`. Assumes forward-forward orientation
/// (bracket placement `t]p]` for record 1 and `]p]t` for record 2), which is
/// the correct representation for a 5′→3′ fusion where partner A contributes
/// the upstream sequence and partner B the downstream sequence.
///
/// The `writer` target must already contain a VCF header (use `write_vcf` for
/// full-pipeline output, which writes the header automatically).
///
/// # Example
/// ```
/// use kam_call::output::write_fusion_bnd_records;
/// use kam_call::fusion::{FusionCall, GenomicLocus};
/// use kam_call::caller::VariantFilter;
///
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
/// let mut buf = Vec::new();
/// write_fusion_bnd_records(&call, &mut buf).unwrap();
/// let text = String::from_utf8(buf).unwrap();
/// assert!(text.contains("chr22"));
/// assert!(text.contains("chr9"));
/// assert!(text.contains("SVTYPE=BND"));
/// assert!(text.contains("MATEID="));
/// ```
pub fn write_fusion_bnd_records(call: &FusionCall, writer: &mut dyn Write) -> io::Result<()> {
    let filter_str = filter_to_str(call.filter);
    let id_prefix = format!("bnd_{}", call.name);
    let id_1 = format!("{id_prefix}_1");
    let id_2 = format!("{id_prefix}_2");

    // Use a placeholder anchor base N when no sequence is available.
    let ref_a = "N";
    let ref_b = "N";

    let base_info = format!(
        "VAF={:.6};VAF_LO={:.6};VAF_HI={:.6};NALT={};NDUPALT={};CONF={:.6};SVTYPE=BND",
        call.vaf,
        call.vaf_ci_low,
        call.vaf_ci_high,
        call.n_molecules,
        call.n_duplex,
        call.confidence,
    );

    // Record 1: partner A breakpoint. ALT = ref_base]chrom_b:pos_b]
    // (forward-forward orientation: closing bracket on right).
    writeln!(
        writer,
        "{chrom_a}\t{pos_a}\t{id1}\t{ref_a}\t{ref_a}]{chrom_b}:{pos_b}]\t.\t{filter}\t{info};MATEID={id2}",
        chrom_a = call.locus_a.chrom,
        pos_a = call.locus_a.start + 1, // convert 0-based to 1-based VCF POS
        id1 = id_1,
        ref_a = ref_a,
        chrom_b = call.locus_b.chrom,
        pos_b = call.locus_b.start + 1,
        filter = filter_str,
        info = base_info,
        id2 = id_2,
    )?;

    // Record 2: partner B breakpoint. ALT = ]chrom_a:pos_a]ref_base
    // (forward-forward orientation: closing bracket on left).
    writeln!(
        writer,
        "{chrom_b}\t{pos_b}\t{id2}\t{ref_b}\t]{chrom_a}:{pos_a}]{ref_b}\t.\t{filter}\t{info};MATEID={id1}",
        chrom_b = call.locus_b.chrom,
        pos_b = call.locus_b.start + 1,
        id2 = id_2,
        ref_b = ref_b,
        chrom_a = call.locus_a.chrom,
        pos_a = call.locus_a.start + 1,
        filter = filter_str,
        info = base_info,
        id1 = id_1,
    )?;

    Ok(())
}

// ─── Private helpers ──────────────────────────────────────────────────────────

/// Format an `Option<T>` as its string value or `"."` when absent.
fn opt_to_str<T: std::fmt::Display>(v: Option<T>) -> String {
    match v {
        Some(x) => x.to_string(),
        None => ".".to_string(),
    }
}

/// Write a delimited (TSV or CSV) format.
///
/// `ml_filter_threshold` is the probability at or above which a call is
/// labelled `ML_PASS`. Callers that do not use a model should pass `0.5`.
fn write_delimited(
    calls: &[VariantCall],
    sep: char,
    writer: &mut dyn Write,
    ml_filter_threshold: f64,
) -> io::Result<()> {
    let s = sep;
    writeln!(
        writer,
        "target_id{s}variant_type{s}ref_seq{s}alt_seq{s}vaf{s}vaf_ci_low{s}vaf_ci_high{s}\
         n_molecules_ref{s}n_molecules_alt{s}n_duplex_alt{s}n_simplex_alt{s}\
         strand_bias_p{s}confidence{s}filter{s}\
         n_simplex_fwd_alt{s}n_simplex_rev_alt{s}n_duplex_ref{s}n_simplex_ref{s}\
         mean_alt_error_prob{s}min_variant_specific_duplex{s}mean_variant_specific_molecules{s}\
         ml_prob{s}ml_filter{s}\
         call_source{s}rescue_min_alt_molecules{s}rescue_alt_duplex{s}rescue_approx_vaf{s}rescue_kmers_found{s}rescue_kmers_total"
    )?;
    for call in calls {
        let ml_prob_str = match call.ml_prob {
            Some(p) => format!("{:.4}", p),
            None => ".".to_string(),
        };
        let ml_filter_str = match call.ml_prob {
            Some(p) if f64::from(p) >= ml_filter_threshold => "ML_PASS",
            Some(_) => "ML_FILTER",
            None => ".",
        };
        let rescue_min = opt_to_str(call.rescue_min_alt_molecules);
        let rescue_dup = opt_to_str(call.rescue_alt_duplex);
        let rescue_vaf = match call.rescue_approx_vaf {
            Some(v) => format!("{:.6}", v),
            None => ".".to_string(),
        };
        let rescue_kf = opt_to_str(call.rescue_kmers_found);
        let rescue_kt = opt_to_str(call.rescue_kmers_total);
        writeln!(
            writer,
            "{}{s}{}{s}{}{s}{}{s}{:.6}{s}{:.6}{s}{:.6}{s}{}{s}{}{s}{}{s}{}{s}{:.6}{s}{:.6}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{:.6}{s}{}{s}{:.4}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}",
            call.target_id,
            variant_type_to_str(call.variant_type),
            bytes_to_str(&call.ref_sequence),
            bytes_to_str(&call.alt_sequence),
            call.vaf,
            call.vaf_ci_low,
            call.vaf_ci_high,
            call.n_molecules_ref,
            call.n_molecules_alt,
            call.n_duplex_alt,
            call.n_simplex_alt,
            call.strand_bias_p,
            call.confidence,
            filter_to_str(call.filter),
            call.n_simplex_fwd_alt,
            call.n_simplex_rev_alt,
            call.n_duplex_ref,
            call.n_simplex_ref,
            call.mean_alt_error_prob,
            call.min_variant_specific_duplex,
            call.mean_variant_specific_molecules,
            ml_prob_str,
            ml_filter_str,
            call.call_source,
            rescue_min,
            rescue_dup,
            rescue_vaf,
            rescue_kf,
            rescue_kt,
        )?;
    }
    Ok(())
}

/// Write calls as a JSON array using `ml_filter_threshold` for the `ml_filter` field.
fn write_json_with_threshold(
    calls: &[VariantCall],
    writer: &mut dyn Write,
    ml_filter_threshold: f64,
) -> io::Result<()> {
    let items: Vec<Value> = calls
        .iter()
        .map(|c| call_to_json(c, ml_filter_threshold))
        .collect();
    let serialised = serde_json::to_string_pretty(&items).map_err(io::Error::other)?;
    writeln!(writer, "{serialised}")
}

/// Serialise a single [`VariantCall`] to a JSON [`Value`].
fn call_to_json(call: &VariantCall, ml_filter_threshold: f64) -> Value {
    let ml_filter = match call.ml_prob {
        Some(p) if f64::from(p) >= ml_filter_threshold => "ML_PASS",
        Some(_) => "ML_FILTER",
        None => ".",
    };
    json!({
        "target_id": call.target_id,
        "variant_type": variant_type_to_str(call.variant_type),
        "ref_seq": bytes_to_str(&call.ref_sequence),
        "alt_seq": bytes_to_str(&call.alt_sequence),
        "vaf": call.vaf,
        "vaf_ci_low": call.vaf_ci_low,
        "vaf_ci_high": call.vaf_ci_high,
        "n_molecules_ref": call.n_molecules_ref,
        "n_molecules_alt": call.n_molecules_alt,
        "n_duplex_alt": call.n_duplex_alt,
        "n_simplex_alt": call.n_simplex_alt,
        "strand_bias_p": call.strand_bias_p,
        "confidence": call.confidence,
        "filter": filter_to_str(call.filter),
        "n_simplex_fwd_alt": call.n_simplex_fwd_alt,
        "n_simplex_rev_alt": call.n_simplex_rev_alt,
        "n_duplex_ref": call.n_duplex_ref,
        "n_simplex_ref": call.n_simplex_ref,
        "mean_alt_error_prob": call.mean_alt_error_prob,
        "min_variant_specific_duplex": call.min_variant_specific_duplex,
        "mean_variant_specific_molecules": call.mean_variant_specific_molecules,
        "ml_prob": call.ml_prob,
        "ml_filter": ml_filter,
        "call_source": call.call_source.to_string(),
        "rescue_min_alt_molecules": call.rescue_min_alt_molecules,
        "rescue_alt_duplex": call.rescue_alt_duplex,
        "rescue_approx_vaf": call.rescue_approx_vaf,
        "rescue_kmers_found": call.rescue_kmers_found,
        "rescue_kmers_total": call.rescue_kmers_total,
    })
}

fn variant_type_to_str(vt: VariantType) -> &'static str {
    match vt {
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
    }
}

fn filter_to_str(vf: VariantFilter) -> &'static str {
    match vf {
        VariantFilter::Pass => "PASS",
        VariantFilter::StrandBias => "StrandBias",
        VariantFilter::LowConfidence => "LowConfidence",
        VariantFilter::LowDuplex => "LowDuplex",
        VariantFilter::CollisionRisk => "CollisionRisk",
        VariantFilter::HighVaf => "HighVaf",
        VariantFilter::NotTargeted => "NotTargeted",
    }
}

/// Convert a byte sequence to a lossy UTF-8 string.
fn bytes_to_str(seq: &[u8]) -> String {
    String::from_utf8_lossy(seq).into_owned()
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::caller::{VariantCall, VariantFilter, VariantType};

    fn make_call(target: &str, ref_seq: &[u8], alt_seq: &[u8]) -> VariantCall {
        use crate::caller::CallSource;
        VariantCall {
            target_id: target.to_owned(),
            variant_type: VariantType::Snv,
            ref_sequence: ref_seq.to_vec(),
            alt_sequence: alt_seq.to_vec(),
            vaf: 0.01,
            vaf_ci_low: 0.005,
            vaf_ci_high: 0.02,
            n_molecules_ref: 990,
            n_molecules_alt: 10,
            n_duplex_alt: 5,
            n_simplex_alt: 5,
            n_simplex_fwd_alt: 3,
            n_simplex_rev_alt: 2,
            n_duplex_ref: 100,
            n_simplex_ref: 890,
            mean_alt_error_prob: 0.001,
            min_variant_specific_duplex: 4,
            mean_variant_specific_molecules: 9.5,
            confidence: 0.999,
            strand_bias_p: 0.8,
            filter: VariantFilter::Pass,
            ml_prob: None,
            call_source: CallSource::Called,
            rescue_min_alt_molecules: None,
            rescue_alt_duplex: None,
            rescue_approx_vaf: None,
            rescue_kmers_found: None,
            rescue_kmers_total: None,
        }
    }

    // Test 1: write_tsv produces correct header and one data line.
    #[test]
    fn tsv_has_header_and_data_line() {
        let calls = vec![make_call("TP53_exon7", b"A", b"T")];
        let mut buf = Vec::new();
        write_tsv(&calls, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 2, "expected header + 1 data line");
        assert!(lines[0].starts_with("target_id\t"));
        assert!(lines[1].contains("TP53_exon7"));
        assert!(lines[1].contains('\t'));
    }

    // Test 2: write_csv uses commas instead of tabs.
    #[test]
    fn csv_uses_commas() {
        let calls = vec![make_call("TP53_exon7", b"A", b"T")];
        let mut buf = Vec::new();
        write_csv(&calls, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains(','), "CSV output must contain commas");
        assert!(!text.contains('\t'), "CSV output must not contain tabs");
    }

    // Test 3: write_json produces valid JSON array.
    #[test]
    fn json_produces_valid_array() {
        let calls = vec![make_call("TP53_exon7", b"A", b"T")];
        let mut buf = Vec::new();
        write_json(&calls, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&text).expect("valid JSON");
        assert!(parsed.is_array());
        assert_eq!(parsed.as_array().unwrap().len(), 1);
        let obj = &parsed[0];
        assert_eq!(obj["target_id"], "TP53_exon7");
        assert_eq!(obj["variant_type"], "SNV");
    }

    // Test 4: write_vcf produces valid VCF header and data line.
    #[test]
    fn vcf_has_header_and_data_line() {
        let calls = vec![make_call("TP53_exon7", b"A", b"T")];
        let mut buf = Vec::new();
        write_vcf(&calls, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("##fileformat=VCFv4.3"));
        assert!(text.contains("#CHROM\tPOS\t"));
        let data_lines: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data_lines.len(), 1);
        assert!(data_lines[0].starts_with("TP53_exon7\t1\t"));
    }

    // Test 5: Multiple variant calls produce multiple lines.
    #[test]
    fn multiple_calls_produce_multiple_lines() {
        let calls = vec![
            make_call("KRAS_G12", b"C", b"A"),
            make_call("TP53_R175", b"G", b"T"),
        ];
        let mut buf = Vec::new();
        write_tsv(&calls, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data_lines: Vec<&str> = text.lines().skip(1).collect();
        assert_eq!(data_lines.len(), 2);
        assert!(data_lines[0].contains("KRAS_G12"));
        assert!(data_lines[1].contains("TP53_R175"));
    }

    // Test 6: Empty calls list → header only (TSV/CSV), empty array (JSON), header only (VCF).
    #[test]
    fn empty_calls_produces_header_only() {
        let calls: Vec<VariantCall> = vec![];

        let mut tsv_buf = Vec::new();
        write_tsv(&calls, &mut tsv_buf).unwrap();
        let tsv = String::from_utf8(tsv_buf).unwrap();
        assert_eq!(tsv.lines().count(), 1, "TSV: only header expected");

        let mut csv_buf = Vec::new();
        write_csv(&calls, &mut csv_buf).unwrap();
        let csv = String::from_utf8(csv_buf).unwrap();
        assert_eq!(csv.lines().count(), 1, "CSV: only header expected");

        let mut json_buf = Vec::new();
        write_json(&calls, &mut json_buf).unwrap();
        let json_text = String::from_utf8(json_buf).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&json_text).unwrap();
        assert_eq!(
            parsed.as_array().unwrap().len(),
            0,
            "JSON: empty array expected"
        );

        let mut vcf_buf = Vec::new();
        write_vcf(&calls, &mut vcf_buf).unwrap();
        let vcf = String::from_utf8(vcf_buf).unwrap();
        let data_lines: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data_lines.len(), 0, "VCF: no data lines expected");
    }

    // Test 7: Special characters in sequences don't break formatting.
    #[test]
    fn special_chars_in_sequences_are_handled() {
        // Sequences with non-ACGT bytes that are still valid UTF-8.
        let call = make_call("target1", b"ACGT", b"ACGN");
        let mut buf = Vec::new();
        write_tsv(&[call], &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("ACGN"));
    }

    // Test 8: VCF uses proper genomic coordinates when target_id is chrN:START-END.
    #[test]
    fn vcf_uses_genomic_coordinates_for_snv() {
        // target_id = "chr2:1000-1100", ref = full 101-bp sequence with C at index 50.
        let mut ref_seq = vec![b'A'; 101];
        let mut alt_seq = ref_seq.clone();
        ref_seq[50] = b'C';
        alt_seq[50] = b'T';
        let mut call = make_call("chr2:1000-1100", &ref_seq, &alt_seq);
        call.variant_type = VariantType::Snv;
        let mut buf = Vec::new();
        write_vcf(&[call], &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let data: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data.len(), 1);
        let fields: Vec<&str> = data[0].split('\t').collect();
        assert_eq!(fields[0], "chr2", "CHROM");
        assert_eq!(fields[1], "1050", "POS = start(1000) + index(50)");
        assert_eq!(fields[3], "C", "REF");
        assert_eq!(fields[4], "T", "ALT");
    }

    // Test 9: write_vcf uses symbolic allele for LargeDeletion.
    #[test]
    fn vcf_large_deletion_uses_symbolic_allele() {
        // Simulate a deletion junction target: ref is longer than alt by 60 bp.
        let ref_seq: Vec<u8> = b"ACGT".repeat(25).to_vec(); // 100 bp
        let alt_seq: Vec<u8> = b"ACGT".repeat(10).to_vec(); // 40 bp
        let mut call = make_call("chr5:1000-1100", &ref_seq, &alt_seq);
        call.variant_type = VariantType::LargeDeletion;
        let mut buf = Vec::new();
        write_vcf(&[call], &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        // Must contain the symbolic allele.
        assert!(text.contains("<DEL>"), "expected <DEL> in VCF output");
        // Must contain SVTYPE and SVLEN in INFO.
        assert!(text.contains("SVTYPE=DEL"), "expected SVTYPE=DEL");
        assert!(text.contains("SVLEN="), "expected SVLEN");
        // Must contain the ALT header definition.
        assert!(text.contains("##ALT=<ID=DEL"), "expected ##ALT=<ID=DEL");
    }

    // ─── SV-EXP-010: VCF/TSV output for new SV types ─────────────────────────

    fn make_sv_call(vt: VariantType, ref_seq: &[u8], alt_seq: &[u8]) -> VariantCall {
        let mut call = make_call("chrX:1000-2000", ref_seq, alt_seq);
        call.variant_type = vt;
        call
    }

    // Test 10a: VCF output for InvDel has SVTYPE=INVDEL and symbolic allele <INVDEL>.
    #[test]
    fn vcf_invdel_has_correct_svtype() {
        let ref_seq: Vec<u8> = b"ACGT".repeat(40).to_vec(); // 160 bp
        let alt_seq: Vec<u8> = b"ACGT".repeat(27).to_vec(); // 108 bp — 52 bp shorter
        let call = make_sv_call(VariantType::InvDel, &ref_seq, &alt_seq);
        let mut buf = Vec::new();
        write_vcf(&[call], &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("<INVDEL>"), "expected <INVDEL> in VCF ALT");
        assert!(
            text.contains("SVTYPE=INVDEL"),
            "expected SVTYPE=INVDEL in INFO"
        );
        assert!(text.contains("SVLEN="), "expected SVLEN in INFO");
        assert!(
            text.contains("##ALT=<ID=INVDEL"),
            "expected INVDEL ALT header"
        );
    }

    // Test 10b: VCF output for NovelInsertion has SVTYPE=INS and symbolic allele <INS>.
    #[test]
    fn vcf_novel_insertion_has_correct_svtype() {
        let ref_seq: Vec<u8> = b"ACGT".repeat(25).to_vec(); // 100 bp
                                                            // alt is longer by 60 bp — a novel insertion.
        let mut alt_seq = ref_seq[..50].to_vec();
        alt_seq.extend_from_slice(&[b'C'; 60]);
        alt_seq.extend_from_slice(&ref_seq[50..]);
        let call = make_sv_call(VariantType::NovelInsertion, &ref_seq, &alt_seq);
        let mut buf = Vec::new();
        write_vcf(&[call], &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("<INS>"), "expected <INS> in VCF ALT");
        assert!(text.contains("SVTYPE=INS"), "expected SVTYPE=INS in INFO");
        assert!(text.contains("##ALT=<ID=INS"), "expected INS ALT header");
    }

    // Test 10c: VCF output for Fusion has SVTYPE=BND and symbolic allele <BND>.
    #[test]
    fn vcf_fusion_has_correct_svtype() {
        let ref_seq: Vec<u8> = b"NNNN".repeat(25).to_vec(); // 100 bp placeholder
        let alt_seq: Vec<u8> = b"NNNN".repeat(25).to_vec(); // same length placeholder
        let call = make_sv_call(VariantType::Fusion, &ref_seq, &alt_seq);
        let mut buf = Vec::new();
        write_vcf(&[call], &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("<BND>"), "expected <BND> in VCF ALT");
        assert!(text.contains("SVTYPE=BND"), "expected SVTYPE=BND in INFO");
        assert!(text.contains("##ALT=<ID=BND"), "expected BND ALT header");
    }

    // Test 10d: TSV output has correct variant_type strings for new SV types.
    #[test]
    fn tsv_new_sv_types_have_correct_type_strings() {
        let ref_seq: Vec<u8> = b"ACGT".repeat(40).to_vec();
        let alt_seq_del: Vec<u8> = b"ACGT".repeat(27).to_vec();

        let calls = vec![
            make_sv_call(VariantType::InvDel, &ref_seq, &alt_seq_del),
            make_sv_call(VariantType::NovelInsertion, &ref_seq, &ref_seq), // same-length placeholder
            make_sv_call(VariantType::Fusion, &ref_seq, &ref_seq),
        ];

        let mut buf = Vec::new();
        write_tsv(&calls, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();

        assert!(
            text.contains("InvDel"),
            "TSV must contain 'InvDel' for InvDel variant"
        );
        assert!(
            text.contains("NovelInsertion"),
            "TSV must contain 'NovelInsertion' for NovelInsertion variant"
        );
        assert!(
            text.contains("Fusion"),
            "TSV must contain 'Fusion' for Fusion variant"
        );
    }

    // Test 10: write_variants dispatches to the correct format.
    #[test]
    fn write_variants_dispatches_correctly() {
        let calls = vec![make_call("test", b"A", b"T")];

        let mut tsv_buf = Vec::new();
        write_variants(&calls, OutputFormat::Tsv, &mut tsv_buf, 0.5).unwrap();
        assert!(String::from_utf8(tsv_buf).unwrap().contains('\t'));

        let mut csv_buf = Vec::new();
        write_variants(&calls, OutputFormat::Csv, &mut csv_buf, 0.5).unwrap();
        assert!(String::from_utf8(csv_buf).unwrap().contains(','));

        let mut json_buf = Vec::new();
        write_variants(&calls, OutputFormat::Json, &mut json_buf, 0.5).unwrap();
        let parsed: serde_json::Value =
            serde_json::from_str(&String::from_utf8(json_buf).unwrap()).unwrap();
        assert!(parsed.is_array());

        let mut vcf_buf = Vec::new();
        write_variants(&calls, OutputFormat::Vcf, &mut vcf_buf, 0.5).unwrap();
        assert!(String::from_utf8(vcf_buf).unwrap().contains("##fileformat"));
    }

    // Test 11: ml_filter column respects the supplied threshold, not the hardcoded 0.5.
    #[test]
    fn ml_filter_respects_custom_threshold() {
        let mut call = make_call("TARGET1", b"A", b"T");
        call.ml_prob = Some(0.46_f32);

        let mut buf = Vec::new();
        write_delimited(&[call.clone()], '\t', &mut buf, 0.449).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(
            output.contains("ML_PASS"),
            "0.46 >= 0.449 should be ML_PASS, got: {output}"
        );

        let mut buf2 = Vec::new();
        write_delimited(&[call], '\t', &mut buf2, 0.5).unwrap();
        let output2 = String::from_utf8(buf2).unwrap();
        assert!(
            output2.contains("ML_FILTER"),
            "0.46 < 0.5 should be ML_FILTER, got: {output2}"
        );
    }

    // ── Additional edge-case tests ───────────────────────────────────────────

    // Test 11: TSV header contains all expected column names in the correct order.
    // Column order matters because downstream parsers (e.g. Python pandas) rely
    // on positional indexing.
    #[test]
    fn tsv_header_has_all_columns_in_order() {
        let mut buf = Vec::new();
        write_tsv(&[], &mut buf).expect("write should succeed");
        let text = String::from_utf8(buf).expect("valid UTF-8");
        let header = text.lines().next().expect("at least one line");
        let cols: Vec<&str> = header.split('\t').collect();
        assert_eq!(cols[0], "target_id", "first column");
        assert_eq!(cols[1], "variant_type", "second column");
        assert_eq!(cols[2], "ref_seq", "third column");
        assert_eq!(cols[3], "alt_seq", "fourth column");
        assert_eq!(cols[4], "vaf", "fifth column");
        // Check key later columns by name presence rather than exact index.
        assert!(
            cols.contains(&"ml_prob"),
            "header must contain ml_prob column"
        );
        assert!(
            cols.contains(&"ml_filter"),
            "header must contain ml_filter column"
        );
        assert!(
            cols.contains(&"call_source"),
            "header must contain call_source column"
        );
        assert!(
            cols.contains(&"rescue_min_alt_molecules"),
            "header must contain rescue column"
        );
    }

    // Test 12: ml_filter column says "." when ml_prob is None (no model loaded).
    // This is the default case when the ML feature is not used.
    #[test]
    fn tsv_ml_filter_dot_when_no_model() {
        let call = make_call("t1", b"A", b"T");
        assert!(call.ml_prob.is_none(), "precondition: no ml_prob");
        let mut buf = Vec::new();
        write_tsv(&[call], &mut buf).expect("write");
        let text = String::from_utf8(buf).expect("UTF-8");
        let data_line = text.lines().nth(1).expect("data line");
        let fields: Vec<&str> = data_line.split('\t').collect();
        // ml_prob and ml_filter are columns 21 and 22 (0-indexed).
        assert_eq!(fields[21], ".", "ml_prob should be '.' when None");
        assert_eq!(fields[22], ".", "ml_filter should be '.' when None");
    }

    // Test 13: ml_filter says ML_PASS when ml_prob >= 0.5 (default threshold).
    // The hardcoded threshold in write_delimited is 0.5.
    #[test]
    fn tsv_ml_filter_pass_at_default_threshold() {
        let mut call = make_call("t1", b"A", b"T");
        call.ml_prob = Some(0.5);
        let mut buf = Vec::new();
        write_tsv(&[call], &mut buf).expect("write");
        let text = String::from_utf8(buf).expect("UTF-8");
        let data_line = text.lines().nth(1).expect("data line");
        let fields: Vec<&str> = data_line.split('\t').collect();
        assert_eq!(
            fields[22], "ML_PASS",
            "ml_prob=0.5 must produce ML_PASS at default threshold"
        );
    }

    // Test 14: ml_filter says ML_FILTER when ml_prob < 0.5 (default threshold).
    #[test]
    fn tsv_ml_filter_fail_below_threshold() {
        let mut call = make_call("t1", b"A", b"T");
        call.ml_prob = Some(0.499);
        let mut buf = Vec::new();
        write_tsv(&[call], &mut buf).expect("write");
        let text = String::from_utf8(buf).expect("UTF-8");
        let data_line = text.lines().nth(1).expect("data line");
        let fields: Vec<&str> = data_line.split('\t').collect();
        assert_eq!(
            fields[22], "ML_FILTER",
            "ml_prob=0.499 must produce ML_FILTER at default threshold"
        );
    }

    // Test 15: JSON output includes all expected keys and respects ml_filter.
    // Verifying the JSON schema matters because external tools parse it.
    #[test]
    fn json_output_has_all_keys_and_ml_filter() {
        let mut call = make_call("gene1", b"C", b"A");
        call.ml_prob = Some(0.8);
        let mut buf = Vec::new();
        write_json(&[call], &mut buf).expect("write");
        let text = String::from_utf8(buf).expect("UTF-8");
        let parsed: serde_json::Value = serde_json::from_str(&text).expect("valid JSON");
        let obj = &parsed[0];
        // Check essential keys exist.
        assert!(obj.get("target_id").is_some(), "must have target_id");
        assert!(obj.get("vaf").is_some(), "must have vaf");
        assert!(obj.get("confidence").is_some(), "must have confidence");
        assert!(obj.get("ml_prob").is_some(), "must have ml_prob");
        assert!(obj.get("ml_filter").is_some(), "must have ml_filter");
        assert!(obj.get("call_source").is_some(), "must have call_source");
        // Verify ml_filter value.
        assert_eq!(
            obj["ml_filter"], "ML_PASS",
            "ml_prob=0.8 should give ML_PASS"
        );
    }

    // Test 16: VCF header contains required lines: ##fileformat, ##INFO, #CHROM.
    #[test]
    fn vcf_header_required_lines() {
        let mut buf = Vec::new();
        write_vcf(&[], &mut buf).expect("write");
        let text = String::from_utf8(buf).expect("UTF-8");
        assert!(
            text.contains("##fileformat=VCFv4.3"),
            "must contain fileformat"
        );
        assert!(text.contains("##INFO=<ID=VAF"), "must contain VAF INFO");
        assert!(
            text.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"),
            "must contain column header"
        );
        assert!(
            text.contains("##FILTER=<ID=PASS"),
            "must contain PASS filter definition"
        );
    }

    // Test 17: VCF FILTER column says "PASS" for passing calls.
    #[test]
    fn vcf_filter_pass_for_passing_call() {
        let call = make_call("simple_target", b"A", b"T");
        assert_eq!(call.filter, VariantFilter::Pass, "precondition");
        let mut buf = Vec::new();
        write_vcf(&[call], &mut buf).expect("write");
        let text = String::from_utf8(buf).expect("UTF-8");
        let data: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data.len(), 1, "one data line");
        let fields: Vec<&str> = data[0].split('\t').collect();
        assert_eq!(fields[6], "PASS", "FILTER column must say PASS");
    }

    // Test 18: VCF FILTER column shows filter name for filtered calls.
    #[test]
    fn vcf_filter_shows_filter_name() {
        let mut call = make_call("simple_target", b"A", b"T");
        call.filter = VariantFilter::StrandBias;
        let mut buf = Vec::new();
        write_vcf(&[call], &mut buf).expect("write");
        let text = String::from_utf8(buf).expect("UTF-8");
        let data: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        let fields: Vec<&str> = data[0].split('\t').collect();
        assert_eq!(
            fields[6], "StrandBias",
            "FILTER column must show StrandBias"
        );
    }

    // Test 19: CSV uses commas (not tabs) and has the same columns as TSV.
    // This is a structural test: CSV and TSV must have identical schemas.
    #[test]
    fn csv_same_column_count_as_tsv() {
        let calls = vec![make_call("t1", b"A", b"T")];
        let mut tsv_buf = Vec::new();
        write_tsv(&calls, &mut tsv_buf).expect("tsv");
        let tsv_text = String::from_utf8(tsv_buf).expect("UTF-8");
        let tsv_cols = tsv_text
            .lines()
            .next()
            .expect("header")
            .split('\t')
            .count();

        let mut csv_buf = Vec::new();
        write_csv(&calls, &mut csv_buf).expect("csv");
        let csv_text = String::from_utf8(csv_buf).expect("UTF-8");
        let csv_cols = csv_text
            .lines()
            .next()
            .expect("header")
            .split(',')
            .count();

        assert_eq!(
            tsv_cols, csv_cols,
            "TSV and CSV must have the same number of columns"
        );
    }

    // Test 20: Call with all optional fields None (rescue fields) does not panic
    // and fills those columns with ".".
    #[test]
    fn call_with_all_optional_none_outputs_dots() {
        let call = make_call("t1", b"A", b"T");
        // make_call already sets rescue_* to None.
        assert!(call.rescue_min_alt_molecules.is_none(), "precondition");
        let mut buf = Vec::new();
        write_tsv(&[call], &mut buf).expect("write should not panic");
        let text = String::from_utf8(buf).expect("UTF-8");
        let data_line = text.lines().nth(1).expect("data line");
        // rescue columns are the last 5 fields.
        let fields: Vec<&str> = data_line.split('\t').collect();
        let n = fields.len();
        for (offset, field) in fields[n - 5..].iter().enumerate() {
            assert_eq!(
                *field, ".",
                "rescue field at offset {offset} from end must be '.' when None"
            );
        }
    }

    // Test 21: SV call (LargeDeletion) produces correct SVLEN in VCF INFO.
    // SVLEN must be negative for deletions per VCF spec.
    #[test]
    fn vcf_large_deletion_svlen_is_negative() {
        let ref_seq: Vec<u8> = vec![b'A'; 100];
        let alt_seq: Vec<u8> = vec![b'A'; 40]; // 60 bp deletion
        let mut call = make_call("chrX:1000-1099", &ref_seq, &alt_seq);
        call.variant_type = VariantType::LargeDeletion;
        let mut buf = Vec::new();
        write_vcf(&[call], &mut buf).expect("write");
        let text = String::from_utf8(buf).expect("UTF-8");
        // SVLEN should be -60.
        assert!(
            text.contains("SVLEN=-60"),
            "LargeDeletion SVLEN must be negative: {text}"
        );
    }

    // Test 22: Fusion call produces two BND records with MATEID linking them.
    // This tests the dedicated write_fusion_bnd_records function.
    #[test]
    fn fusion_bnd_records_linked_by_mateid() {
        let fusion_call = FusionCall {
            name: "BCR_ABL1".to_string(),
            locus_a: crate::fusion::GenomicLocus {
                chrom: "chr22".to_string(),
                start: 23_632_500,
                end: 23_632_550,
            },
            locus_b: crate::fusion::GenomicLocus {
                chrom: "chr9".to_string(),
                start: 130_854_000,
                end: 130_854_050,
            },
            vaf: 0.01,
            vaf_ci_low: 0.005,
            vaf_ci_high: 0.02,
            n_molecules: 10,
            n_duplex: 5,
            confidence: 0.999,
            filter: VariantFilter::Pass,
        };
        let mut buf = Vec::new();
        write_fusion_bnd_records(&fusion_call, &mut buf).expect("write");
        let text = String::from_utf8(buf).expect("UTF-8");
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 2, "must produce exactly two BND records");
        // First record: chr22, second: chr9.
        assert!(lines[0].starts_with("chr22\t"), "first record is partner A");
        assert!(lines[1].starts_with("chr9\t"), "second record is partner B");
        // Both records have SVTYPE=BND.
        assert!(lines[0].contains("SVTYPE=BND"), "record 1 SVTYPE");
        assert!(lines[1].contains("SVTYPE=BND"), "record 2 SVTYPE");
        // MATEIDs cross-reference each other.
        assert!(
            lines[0].contains("MATEID=bnd_BCR_ABL1_2"),
            "record 1 MATEID must point to record 2"
        );
        assert!(
            lines[1].contains("MATEID=bnd_BCR_ABL1_1"),
            "record 2 MATEID must point to record 1"
        );
    }

    // Test 23: JSON output for an empty call list is exactly "[]".
    // Downstream JSON parsers must handle this without error.
    #[test]
    fn json_empty_is_empty_array() {
        let mut buf = Vec::new();
        write_json(&[], &mut buf).expect("write");
        let text = String::from_utf8(buf).expect("UTF-8");
        assert_eq!(text.trim(), "[]", "empty calls must produce '[]'");
    }
}
