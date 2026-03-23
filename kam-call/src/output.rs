//! Writers for variant calls in TSV, CSV, JSON, and VCF formats.
//!
//! All writers accept a `&[VariantCall]` slice and a `&mut dyn Write` target,
//! enabling use with files, sockets, or in-memory buffers equally.

use std::io::{self, Write};

use serde_json::{json, Value};

use crate::allele::extract_minimal_allele;
use crate::caller::{VariantCall, VariantFilter, VariantType};

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
/// # Example
/// ```
/// use kam_call::output::{write_variants, OutputFormat};
/// use kam_call::caller::{VariantCall, VariantType, VariantFilter};
///
/// let calls: Vec<VariantCall> = vec![];
/// let mut buf = Vec::new();
/// write_variants(&calls, OutputFormat::Tsv, &mut buf).unwrap();
/// assert!(!buf.is_empty());
/// ```
pub fn write_variants(
    calls: &[VariantCall],
    format: OutputFormat,
    writer: &mut dyn Write,
) -> io::Result<()> {
    match format {
        OutputFormat::Tsv => write_tsv(calls, writer),
        OutputFormat::Csv => write_csv(calls, writer),
        OutputFormat::Json => write_json(calls, writer),
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
    write_delimited(calls, '\t', writer)
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
    write_delimited(calls, ',', writer)
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
    let items: Vec<Value> = calls.iter().map(call_to_json).collect();
    let serialised = serde_json::to_string_pretty(&items).map_err(io::Error::other)?;
    writeln!(writer, "{serialised}")
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
        VariantType::LargeDeletion | VariantType::TandemDuplication | VariantType::Inversion
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
        _ => unreachable!("write_sv_vcf_record called with non-SV type"),
    };

    // SVLEN: negative for deletions, positive for duplications, 0 for inversions.
    let svlen: i64 = match call.variant_type {
        VariantType::LargeDeletion => {
            -((call.ref_sequence.len() as i64) - (call.alt_sequence.len() as i64))
        }
        VariantType::TandemDuplication => {
            (call.alt_sequence.len() as i64) - (call.ref_sequence.len() as i64)
        }
        VariantType::Inversion => 0,
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

// ─── Private helpers ──────────────────────────────────────────────────────────

/// Write a delimited (TSV or CSV) format.
fn write_delimited(calls: &[VariantCall], sep: char, writer: &mut dyn Write) -> io::Result<()> {
    let s = sep;
    writeln!(
        writer,
        "target_id{s}variant_type{s}ref_seq{s}alt_seq{s}vaf{s}vaf_ci_low{s}vaf_ci_high{s}\
         n_molecules_ref{s}n_molecules_alt{s}n_duplex_alt{s}n_simplex_alt{s}\
         strand_bias_p{s}confidence{s}filter"
    )?;
    for call in calls {
        writeln!(
            writer,
            "{}{s}{}{s}{}{s}{}{s}{:.6}{s}{:.6}{s}{:.6}{s}{}{s}{}{s}{}{s}{}{s}{:.6}{s}{:.6}{s}{}",
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
        )?;
    }
    Ok(())
}

/// Serialise a single [`VariantCall`] to a JSON [`Value`].
fn call_to_json(call: &VariantCall) -> Value {
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
            confidence: 0.999,
            strand_bias_p: 0.8,
            filter: VariantFilter::Pass,
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

    // Test 10: write_variants dispatches to the correct format.
    #[test]
    fn write_variants_dispatches_correctly() {
        let calls = vec![make_call("test", b"A", b"T")];

        let mut tsv_buf = Vec::new();
        write_variants(&calls, OutputFormat::Tsv, &mut tsv_buf).unwrap();
        assert!(String::from_utf8(tsv_buf).unwrap().contains('\t'));

        let mut csv_buf = Vec::new();
        write_variants(&calls, OutputFormat::Csv, &mut csv_buf).unwrap();
        assert!(String::from_utf8(csv_buf).unwrap().contains(','));

        let mut json_buf = Vec::new();
        write_variants(&calls, OutputFormat::Json, &mut json_buf).unwrap();
        let parsed: serde_json::Value =
            serde_json::from_str(&String::from_utf8(json_buf).unwrap()).unwrap();
        assert!(parsed.is_array());

        let mut vcf_buf = Vec::new();
        write_variants(&calls, OutputFormat::Vcf, &mut vcf_buf).unwrap();
        assert!(String::from_utf8(vcf_buf).unwrap().contains("##fileformat"));
    }
}
