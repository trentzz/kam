//! K-mer rescue probe for TI targets that produce no matching call.
//!
//! For each target variant in the TI set that has no PASS or sub-threshold
//! call, this module queries the k-mer index for alt-supporting evidence.

use std::collections::HashSet;

use kam_core::kmer::KmerIndex;
use kam_index::encode::{canonical, KmerIterator};
use kam_index::HashKmerIndex;

/// Evidence extracted from the k-mer index at a TI target position.
#[derive(Debug, Clone)]
pub struct RescueEvidence {
    /// Minimum molecule count across alt-specific k-mers (lower bound on alt support).
    pub min_alt_molecules: u32,
    /// Sum of duplex counts across alt-specific k-mers (overcount; indicative only).
    pub alt_duplex: u32,
    /// Number of alt-specific k-mers found in the index with >0 molecules.
    pub n_alt_kmers_found: usize,
    /// Total alt-specific k-mers checked.
    pub n_alt_kmers_total: usize,
    /// Mean molecule count across ref-specific k-mers (reference depth estimate).
    pub mean_ref_molecules: f32,
    /// Approximate VAF: min_alt / (min_alt + mean_ref).
    pub approx_vaf: f32,
}

/// Query the k-mer index for alt-specific evidence.
///
/// `ref_seq` and `alt_seq` are full target-window sequences. Returns `None`
/// if either sequence is shorter than `k` or no alt-specific k-mers exist.
///
/// # Example
///
/// ```
/// use kam::rescue::probe_ti_target;
/// ```
pub fn probe_ti_target(
    index: &HashKmerIndex,
    ref_seq: &[u8],
    alt_seq: &[u8],
    k: usize,
) -> Option<RescueEvidence> {
    if ref_seq.len() < k || alt_seq.len() < k {
        return None;
    }

    let ref_kmers: HashSet<u64> = KmerIterator::new(ref_seq, k)
        .map(|(_, kmer)| canonical(kmer, k))
        .collect();
    let alt_kmers: HashSet<u64> = KmerIterator::new(alt_seq, k)
        .map(|(_, kmer)| canonical(kmer, k))
        .collect();

    let alt_specific: Vec<u64> = alt_kmers.difference(&ref_kmers).copied().collect();
    let ref_specific: Vec<u64> = ref_kmers.difference(&alt_kmers).copied().collect();

    let n_alt_kmers_total = alt_specific.len();
    if n_alt_kmers_total == 0 {
        return None;
    }

    let mut min_alt_molecules = u32::MAX;
    let mut alt_duplex_sum = 0u32;
    let mut n_alt_kmers_found = 0usize;

    for &kmer in &alt_specific {
        if let Some(ev) = index.get(kmer) {
            if ev.n_molecules > 0 {
                n_alt_kmers_found += 1;
                min_alt_molecules = min_alt_molecules.min(ev.n_molecules);
                alt_duplex_sum = alt_duplex_sum.saturating_add(ev.n_duplex);
            }
        }
    }
    if min_alt_molecules == u32::MAX {
        min_alt_molecules = 0;
    }

    let mean_ref_molecules = if ref_specific.is_empty() {
        0.0
    } else {
        let sum: u64 = ref_specific
            .iter()
            .filter_map(|&km| index.get(km))
            .map(|ev| ev.n_molecules as u64)
            .sum();
        sum as f32 / ref_specific.len() as f32
    };

    let denom = min_alt_molecules as f32 + mean_ref_molecules;
    let approx_vaf = if denom > 0.0 {
        min_alt_molecules as f32 / denom
    } else {
        0.0
    };

    Some(RescueEvidence {
        min_alt_molecules,
        alt_duplex: alt_duplex_sum,
        n_alt_kmers_found,
        n_alt_kmers_total,
        mean_ref_molecules,
        approx_vaf,
    })
}

/// An alt-sequence entry parsed from the `--alt-as-ref` FASTA.
#[derive(Debug, Clone)]
pub struct AltWalkEntry {
    /// Target window identifier (e.g. `"chr1:1000-1100"`).
    pub target_id: String,
    /// REF allele from the FASTA header `REF=` field.
    pub ref_allele: String,
    /// ALT allele from the FASTA header `ALT=` field.
    pub alt_allele: String,
    /// Pre-built alt-allele window sequence.
    pub alt_seq: Vec<u8>,
}

/// Parse an alt-walk FASTA file into a list of [`AltWalkEntry`] values.
///
/// Expects headers of the form `chrN:start-end REF=xxx ALT=yyy` (the format
/// produced by `generate_alt_seqs_from_norm.py`). Entries whose headers do
/// not contain both `REF=` and `ALT=` tokens are skipped with a warning.
///
/// # Example
///
/// ```
/// use std::io::Write;
/// use kam::rescue::load_alt_walk_entries;
/// let fasta = b">chr1:100-200 REF=G ALT=T\nACGT\n";
/// let tmp = tempfile::NamedTempFile::new().unwrap();
/// tmp.as_file().write_all(fasta).unwrap();
/// let entries = load_alt_walk_entries(tmp.path()).unwrap();
/// assert_eq!(entries.len(), 1);
/// assert_eq!(entries[0].target_id, "chr1:100-200");
/// assert_eq!(entries[0].ref_allele, "G");
/// assert_eq!(entries[0].alt_allele, "T");
/// ```
pub fn load_alt_walk_entries(
    path: &std::path::Path,
) -> Result<Vec<AltWalkEntry>, Box<dyn std::error::Error>> {
    use needletail::parse_fastx_file;
    let mut entries = Vec::new();
    let mut reader = parse_fastx_file(path)?;
    while let Some(result) = reader.next() {
        let rec = result?;
        let header = String::from_utf8_lossy(rec.id());
        // target_id is the part before the first space.
        let mut parts = header.splitn(2, ' ');
        let target_id = parts.next().unwrap_or("").to_string();
        let rest = parts.next().unwrap_or("");
        // Parse REF= and ALT= tokens.
        let ref_allele = rest
            .split_whitespace()
            .find_map(|t| t.strip_prefix("REF="))
            .map(str::to_string);
        let alt_allele = rest
            .split_whitespace()
            .find_map(|t| t.strip_prefix("ALT="))
            .map(str::to_string);
        let (ref_allele, alt_allele) = match (ref_allele, alt_allele) {
            (Some(r), Some(a)) => (r, a),
            _ => {
                eprintln!(
                    "[rescue/alt-walk] skipping entry with unparseable header: {}",
                    header
                );
                continue;
            }
        };
        let alt_seq: Vec<u8> = rec.seq().iter().map(|&b| b.to_ascii_uppercase()).collect();
        entries.push(AltWalkEntry {
            target_id,
            ref_allele,
            alt_allele,
            alt_seq,
        });
    }
    Ok(entries)
}

/// Construct an alt target-window sequence from a ref sequence and a VCF variant.
///
/// - `ref_seq`: full target window reference bytes
/// - `target_start_0based`: 0-based genomic start of the target window
/// - `vcf_pos_1based`: 1-based VCF POS field
/// - `vcf_ref`: VCF REF allele bytes
/// - `vcf_alt`: VCF ALT allele bytes
///
/// Returns `None` if the variant falls outside the window or the ref allele
/// does not match.
///
/// # Example
///
/// ```
/// use kam::rescue::build_alt_seq;
/// let ref_seq = b"ACGTACGT";
/// // SNV at window offset 3 (0-based), target_start=100, vcf_pos=103 (1-based).
/// // multiseqex labels target_start as VCF_POS - flank (effectively 1-based), so
/// // offset = vcf_pos_1based - target_start_0based = 103 - 100 = 3.
/// let alt = build_alt_seq(ref_seq, 100, 103, b"T", b"G");
/// assert_eq!(alt, Some(b"ACGGACGT".to_vec()));
/// ```
pub fn build_alt_seq(
    ref_seq: &[u8],
    target_start_0based: i64,
    vcf_pos_1based: i64,
    vcf_ref: &[u8],
    vcf_alt: &[u8],
) -> Option<Vec<u8>> {
    let offset = vcf_pos_1based - target_start_0based;
    if offset < 0 {
        return None;
    }
    let offset = offset as usize;
    if offset + vcf_ref.len() > ref_seq.len() {
        return None;
    }
    let ref_window = &ref_seq[offset..offset + vcf_ref.len()];
    // Case-insensitive comparison (FASTA can have lowercase masking).
    if !ref_window
        .iter()
        .zip(vcf_ref.iter())
        .all(|(a, b)| a.eq_ignore_ascii_case(b))
    {
        return None;
    }
    let mut alt_seq = Vec::with_capacity(ref_seq.len());
    alt_seq.extend_from_slice(&ref_seq[..offset]);
    alt_seq.extend_from_slice(vcf_alt);
    alt_seq.extend_from_slice(&ref_seq[offset + vcf_ref.len()..]);
    Some(alt_seq)
}

/// Parse a target_id of the form "chrN:start-end" into (chrom, start_0based, end_exclusive).
///
/// Returns `None` if the format is unexpected.
///
/// # Example
///
/// ```
/// use kam::rescue::parse_target_id;
/// let result = parse_target_id("chr1:100-200");
/// assert_eq!(result, Some(("chr1", 100, 200)));
/// ```
pub fn parse_target_id(target_id: &str) -> Option<(&str, i64, i64)> {
    let colon = target_id.find(':')?;
    let chrom = &target_id[..colon];
    let range = &target_id[colon + 1..];
    let dash = range.find('-')?;
    let start: i64 = range[..dash].parse().ok()?;
    let end: i64 = range[dash + 1..].parse().ok()?;
    Some((chrom, start, end))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_target_id_basic() {
        let result = parse_target_id("chr17:7572800-7572900");
        assert_eq!(result, Some(("chr17", 7_572_800, 7_572_900)));
    }

    #[test]
    fn parse_target_id_invalid() {
        assert!(parse_target_id("no_colon").is_none());
        assert!(parse_target_id("chr1:abc-200").is_none());
    }

    #[test]
    fn build_alt_seq_snv() {
        let ref_seq = b"ACGTACGT";
        // target_start=100 (0-based), vcf_pos=103 (1-based) → offset = 103-100 = 3
        // ref_seq[3] = 'T'; substitute with 'G'.
        let alt = build_alt_seq(ref_seq, 100, 103, b"T", b"G");
        assert_eq!(alt, Some(b"ACGGACGT".to_vec()));
    }

    #[test]
    fn build_alt_seq_out_of_window() {
        let ref_seq = b"ACGT";
        // vcf_pos before window start.
        let alt = build_alt_seq(ref_seq, 100, 99, b"A", b"T");
        assert!(alt.is_none());
    }

    #[test]
    fn build_alt_seq_ref_mismatch() {
        let ref_seq = b"ACGT";
        // ref allele does not match window at that position.
        let alt = build_alt_seq(ref_seq, 100, 101, b"T", b"G");
        assert!(alt.is_none());
    }

    #[test]
    fn load_alt_walk_entries_basic() {
        use std::io::Write;
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        writeln!(tmp, ">chr1:100-200 REF=G ALT=T").unwrap();
        writeln!(tmp, "ACGT").unwrap();
        writeln!(tmp, ">chr2:500-600 REF=ACGT ALT=A").unwrap();
        writeln!(tmp, "TTTTACGT").unwrap();
        let entries = load_alt_walk_entries(tmp.path()).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].target_id, "chr1:100-200");
        assert_eq!(entries[0].ref_allele, "G");
        assert_eq!(entries[0].alt_allele, "T");
        assert_eq!(entries[1].ref_allele, "ACGT");
    }
}
