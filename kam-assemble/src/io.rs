//! FASTQ I/O for paired-end Twist UMI duplex reads.
//!
//! Uses `needletail` for streaming FASTQ parsing of R1/R2 file pairs.
//! Both files must contain the same number of records in the same order.

use std::path::Path;

use crate::parser::{parse_read_pair, DropReason, ParseStats, ParsedReadPair, ParserConfig};

// ── Public API ────────────────────────────────────────────────────────────────

/// Read paired FASTQ files and parse all read pairs.
///
/// Opens `r1_path` and `r2_path` with `needletail`, iterates both files in
/// lockstep, and passes each pair through [`crate::parser::parse_read_pair`].
///
/// Reads that pass the parser filters are collected into the returned `Vec`.
/// Reads that are dropped (e.g. too short, low UMI quality) are counted in
/// [`ParseStats`] but not included in the `Vec`.
///
/// # Errors
///
/// Returns an error if:
/// - Either file cannot be opened or contains malformed FASTQ records.
/// - R1 and R2 have a different number of records.
///
/// # Examples
///
/// ```no_run
/// use std::path::Path;
/// use kam_assemble::io::read_fastq_pairs;
/// use kam_assemble::parser::ParserConfig;
///
/// let config = ParserConfig::default();
/// let (pairs, stats) = read_fastq_pairs(
///     Path::new("R1.fastq"),
///     Path::new("R2.fastq"),
///     &config,
/// ).unwrap();
/// println!("Passed: {}", stats.n_passed);
/// ```
pub fn read_fastq_pairs(
    r1_path: &Path,
    r2_path: &Path,
    config: &ParserConfig,
) -> Result<(Vec<ParsedReadPair>, ParseStats), Box<dyn std::error::Error>> {
    // needletail returns an `EmptyFile` error when the file contains no bytes.
    // Treat that as a valid empty input rather than a hard failure.
    let r1_reader_result = needletail::parse_fastx_file(r1_path);
    let r2_reader_result = needletail::parse_fastx_file(r2_path);

    let is_empty_file_err = |e: &needletail::errors::ParseError| {
        e.kind == needletail::errors::ParseErrorKind::EmptyFile
    };

    // If both are empty, return immediately.
    let (r1_empty, r2_empty) = (
        r1_reader_result
            .as_ref()
            .err()
            .map(is_empty_file_err)
            .unwrap_or(false),
        r2_reader_result
            .as_ref()
            .err()
            .map(is_empty_file_err)
            .unwrap_or(false),
    );

    if r1_empty && r2_empty {
        return Ok((Vec::new(), ParseStats::default()));
    }
    if r1_empty || r2_empty {
        return Err("One file is empty and the other is not: mismatched record counts".into());
    }

    let mut r1_reader = r1_reader_result?;
    let mut r2_reader = r2_reader_result?;

    let mut pairs: Vec<ParsedReadPair> = Vec::new();
    let mut stats = ParseStats::default();

    loop {
        let r1_rec = r1_reader.next();
        let r2_rec = r2_reader.next();

        match (r1_rec, r2_rec) {
            // Both files exhausted — normal end.
            (None, None) => break,

            // R1 has a record but R2 does not (R2 is shorter).
            (Some(_), None) => {
                return Err(
                    "R1 has more records than R2: files have mismatched record counts".into(),
                );
            }

            // R2 has a record but R1 does not (R1 is shorter).
            (None, Some(_)) => {
                return Err(
                    "R2 has more records than R1: files have mismatched record counts".into(),
                );
            }

            // Both files have a record — parse the pair.
            (Some(r1_result), Some(r2_result)) => {
                let r1 = r1_result?;
                let r2 = r2_result?;

                let r1_seq = r1.seq();
                let r2_seq = r2.seq();
                // `qual()` returns `Option<&[u8]>`; fall back to an empty
                // slice when quality data is absent (FASTA input).
                let r1_qual: &[u8] = r1.qual().unwrap_or(&[]);
                let r2_qual: &[u8] = r2.qual().unwrap_or(&[]);

                stats.n_processed += 1;

                match parse_read_pair(&r1_seq, r1_qual, &r2_seq, r2_qual, config)? {
                    crate::parser::ParseResult::Ok(parsed) => {
                        stats.n_passed += 1;
                        pairs.push(parsed);
                    }
                    crate::parser::ParseResult::Dropped { reason, .. } => match reason {
                        DropReason::ReadTooShort => stats.n_read_too_short += 1,
                        DropReason::TemplateTooShort => stats.n_template_too_short += 1,
                        DropReason::LowUmiQuality => stats.n_low_umi_quality += 1,
                    },
                }
            }
        }
    }

    Ok((pairs, stats))
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write as IoWrite;

    // ── Helpers ───────────────────────────────────────────────────────────

    /// Write a minimal FASTQ file from a slice of (name, seq, qual) tuples.
    fn write_fastq(path: &Path, records: &[(&str, &str, &str)]) {
        let mut f = fs::File::create(path).expect("create fastq");
        for (name, seq, qual) in records {
            writeln!(f, "@{name}").unwrap();
            writeln!(f, "{seq}").unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{qual}").unwrap();
        }
    }

    /// Build a standard Twist UMI read string: 5 bp UMI + 2 bp skip + template.
    fn make_read_str(umi: &str, skip: &str, template: &str) -> String {
        format!("{umi}{skip}{template}")
    }

    /// Build a quality string of `n` copies of `'I'` (Phred 40).
    fn qual_str(n: usize) -> String {
        "I".repeat(n)
    }

    // ── Test 1: Valid small FASTQ pair produces correct ParsedReadPairs ───

    #[test]
    fn valid_fastq_pair_parses_correctly() {
        let dir = tempdir();
        let r1_path = dir.join("R1.fastq");
        let r2_path = dir.join("R2.fastq");

        let r1_seq = make_read_str("ACGTA", "TG", "NNNNNNNNNNNNNNNNNNNN");
        let r2_seq = make_read_str("TGCAT", "AC", "NNNNNNNNNNNNNNNNNNNN");
        let r1_qual = qual_str(r1_seq.len());
        let r2_qual = qual_str(r2_seq.len());

        write_fastq(&r1_path, &[("read1", &r1_seq, &r1_qual)]);
        write_fastq(&r2_path, &[("read1", &r2_seq, &r2_qual)]);

        let config = ParserConfig::default();
        let (pairs, stats) = read_fastq_pairs(&r1_path, &r2_path, &config).unwrap();

        assert_eq!(pairs.len(), 1, "expected one parsed pair");
        assert_eq!(stats.n_processed, 1);
        assert_eq!(stats.n_passed, 1);
        assert_eq!(stats.n_read_too_short, 0);

        let p = &pairs[0];
        assert_eq!(p.umi_r1, *b"ACGTA");
        assert_eq!(p.umi_r2, *b"TGCAT");
        assert_eq!(p.skip_r1, *b"TG");
        assert_eq!(p.skip_r2, *b"AC");
        assert_eq!(p.template_r1.len(), 20);
        assert_eq!(p.template_r2.len(), 20);
    }

    // ── Test 2: Mismatched record counts returns error ────────────────────

    #[test]
    fn mismatched_record_counts_returns_error() {
        let dir = tempdir();
        let r1_path = dir.join("R1.fastq");
        let r2_path = dir.join("R2.fastq");

        let seq = make_read_str("ACGTA", "TG", "NNNNNNNNNN");
        let qual = qual_str(seq.len());

        // R1 has 2 records, R2 has 1.
        write_fastq(&r1_path, &[("read1", &seq, &qual), ("read2", &seq, &qual)]);
        write_fastq(&r2_path, &[("read1", &seq, &qual)]);

        let config = ParserConfig::default();
        let result = read_fastq_pairs(&r1_path, &r2_path, &config);
        assert!(
            result.is_err(),
            "expected error for mismatched record counts"
        );
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("mismatched"),
            "error message should mention mismatch: {msg}"
        );
    }

    // ── Test 3: Empty files return empty vec with zero stats ──────────────

    #[test]
    fn empty_files_return_empty_vec() {
        let dir = tempdir();
        let r1_path = dir.join("R1.fastq");
        let r2_path = dir.join("R2.fastq");

        write_fastq(&r1_path, &[]);
        write_fastq(&r2_path, &[]);

        let config = ParserConfig::default();
        let (pairs, stats) = read_fastq_pairs(&r1_path, &r2_path, &config).unwrap();

        assert!(pairs.is_empty(), "expected empty vec");
        assert_eq!(stats.n_processed, 0);
        assert_eq!(stats.n_passed, 0);
    }

    // ── Test 4: ParseStats correctly counts drops ─────────────────────────

    #[test]
    fn parse_stats_counts_drops() {
        let dir = tempdir();
        let r1_path = dir.join("R1.fastq");
        let r2_path = dir.join("R2.fastq");

        // Record 1: valid (5 UMI + 2 skip + 20 template = 27 bp)
        let good_seq = make_read_str("ACGTA", "TG", "NNNNNNNNNNNNNNNNNNNN");
        let good_qual = qual_str(good_seq.len());

        // Record 2: too short (only 4 bp — can't extract UMI+skip of 7 bp)
        let short_seq = "ACGT";
        let short_qual = qual_str(short_seq.len());

        write_fastq(
            &r1_path,
            &[
                ("read1", &good_seq, &good_qual),
                ("read2", short_seq, &short_qual),
            ],
        );
        write_fastq(
            &r2_path,
            &[
                ("read1", &good_seq, &good_qual),
                ("read2", short_seq, &short_qual),
            ],
        );

        let config = ParserConfig::default();
        let (pairs, stats) = read_fastq_pairs(&r1_path, &r2_path, &config).unwrap();

        assert_eq!(pairs.len(), 1, "only the valid pair should be collected");
        assert_eq!(stats.n_processed, 2);
        assert_eq!(stats.n_passed, 1);
        assert_eq!(stats.n_read_too_short, 1);
        assert_eq!(stats.n_template_too_short, 0);
        assert_eq!(stats.n_low_umi_quality, 0);
    }

    // ── Test 5: Multiple valid pairs all parsed ───────────────────────────

    #[test]
    fn multiple_valid_pairs_all_collected() {
        let dir = tempdir();
        let r1_path = dir.join("R1.fastq");
        let r2_path = dir.join("R2.fastq");

        let seqs: Vec<(String, String)> = (0..5)
            .map(|i| {
                let seq = make_read_str("ACGTA", "TG", &"N".repeat(10 + i));
                let qual = qual_str(seq.len());
                (seq, qual)
            })
            .collect();

        let r1_records: Vec<(&str, &str, &str)> = seqs
            .iter()
            .enumerate()
            .map(|(i, (s, q))| {
                let name = Box::leak(format!("read{i}").into_boxed_str()) as &str;
                (name, s.as_str(), q.as_str())
            })
            .collect();

        write_fastq(&r1_path, &r1_records);
        write_fastq(&r2_path, &r1_records);

        let config = ParserConfig::default();
        let (pairs, stats) = read_fastq_pairs(&r1_path, &r2_path, &config).unwrap();

        assert_eq!(pairs.len(), 5);
        assert_eq!(stats.n_processed, 5);
        assert_eq!(stats.n_passed, 5);
    }

    // ── Test 6: R2-shorter mismatch also returns error ────────────────────

    #[test]
    fn r2_longer_than_r1_returns_error() {
        let dir = tempdir();
        let r1_path = dir.join("R1.fastq");
        let r2_path = dir.join("R2.fastq");

        let seq = make_read_str("ACGTA", "TG", "NNNNNNNNNN");
        let qual = qual_str(seq.len());

        // R1 has 1 record, R2 has 2.
        write_fastq(&r1_path, &[("read1", &seq, &qual)]);
        write_fastq(&r2_path, &[("read1", &seq, &qual), ("read2", &seq, &qual)]);

        let config = ParserConfig::default();
        let result = read_fastq_pairs(&r1_path, &r2_path, &config);
        assert!(
            result.is_err(),
            "expected error for mismatched record counts"
        );
    }

    // ── Test 7: Gzip-compressed FASTQ pair parses correctly ───────────────

    #[test]
    fn gzip_fastq_pair_parses_correctly() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write as IoWrite2;

        let dir = tempdir();
        let r1_gz_path = dir.join("R1.fastq.gz");
        let r2_gz_path = dir.join("R2.fastq.gz");

        let r1_seq = make_read_str("ACGTA", "TG", "NNNNNNNNNNNNNNNNNNNN");
        let r2_seq = make_read_str("TGCAT", "AC", "NNNNNNNNNNNNNNNNNNNN");
        let r1_qual = qual_str(r1_seq.len());
        let r2_qual = qual_str(r2_seq.len());

        // Write gzip-compressed FASTQ for R1.
        {
            let f = fs::File::create(&r1_gz_path).expect("create R1.fastq.gz");
            let mut gz = GzEncoder::new(f, Compression::default());
            write!(gz, "@read1\n{r1_seq}\n+\n{r1_qual}\n").unwrap();
            gz.finish().expect("finish R1 gzip");
        }

        // Write gzip-compressed FASTQ for R2.
        {
            let f = fs::File::create(&r2_gz_path).expect("create R2.fastq.gz");
            let mut gz = GzEncoder::new(f, Compression::default());
            write!(gz, "@read1\n{r2_seq}\n+\n{r2_qual}\n").unwrap();
            gz.finish().expect("finish R2 gzip");
        }

        let config = ParserConfig::default();
        let (pairs, stats) = read_fastq_pairs(&r1_gz_path, &r2_gz_path, &config).unwrap();

        assert_eq!(pairs.len(), 1, "expected one parsed pair from gzip input");
        assert_eq!(stats.n_processed, 1);
        assert_eq!(stats.n_passed, 1);

        let p = &pairs[0];
        assert_eq!(p.umi_r1, *b"ACGTA");
        assert_eq!(p.umi_r2, *b"TGCAT");
        assert_eq!(p.skip_r1, *b"TG");
        assert_eq!(p.skip_r2, *b"AC");
        assert_eq!(p.template_r1.len(), 20);
        assert_eq!(p.template_r2.len(), 20);
    }

    // ── Helpers ───────────────────────────────────────────────────────────

    /// Create a temporary directory for test files.
    ///
    /// Returns a [`TempDir`] handle — the directory is deleted when it drops.
    fn tempdir() -> TempDir {
        TempDir::new()
    }

    /// A minimal RAII temporary directory.
    struct TempDir {
        path: std::path::PathBuf,
    }

    impl TempDir {
        fn new() -> Self {
            use std::sync::atomic::{AtomicU64, Ordering};
            static COUNTER: AtomicU64 = AtomicU64::new(0);
            let id = COUNTER.fetch_add(1, Ordering::Relaxed);
            let path = std::env::temp_dir()
                .join(format!("kam_assemble_io_test_{id}_{}", std::process::id()));
            fs::create_dir_all(&path).expect("create temp dir");
            Self { path }
        }

        fn join(&self, name: &str) -> std::path::PathBuf {
            self.path.join(name)
        }
    }

    impl Drop for TempDir {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.path);
        }
    }
}
