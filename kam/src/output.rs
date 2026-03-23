//! Shared output-format helpers used by the `call` and `run` subcommands.

use kam_call::output::OutputFormat;

/// Parse a comma-separated format string into [`OutputFormat`] values.
///
/// Unknown tokens produce an error. An empty string defaults to TSV.
///
/// # Errors
///
/// Returns an error string when an unrecognised format token is encountered.
pub fn parse_output_formats(s: &str) -> Result<Vec<OutputFormat>, Box<dyn std::error::Error>> {
    let mut formats = Vec::new();
    for raw in s.split(',') {
        let token = raw.trim().to_ascii_lowercase();
        if token.is_empty() {
            continue;
        }
        let fmt = match token.as_str() {
            "tsv" => OutputFormat::Tsv,
            "csv" => OutputFormat::Csv,
            "json" => OutputFormat::Json,
            "vcf" => OutputFormat::Vcf,
            other => {
                return Err(format!("unknown output format: '{other}'").into());
            }
        };
        formats.push(fmt);
    }
    if formats.is_empty() {
        formats.push(OutputFormat::Tsv);
    }
    Ok(formats)
}

/// Return the file extension string for an output format.
pub fn format_extension(fmt: OutputFormat) -> &'static str {
    match fmt {
        OutputFormat::Tsv => "tsv",
        OutputFormat::Csv => "csv",
        OutputFormat::Json => "json",
        OutputFormat::Vcf => "vcf",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_output_formats_single() {
        let fmts = parse_output_formats("tsv").unwrap();
        assert_eq!(fmts.len(), 1);
        assert_eq!(format_extension(fmts[0]), "tsv");
    }

    #[test]
    fn parse_output_formats_multi() {
        let fmts = parse_output_formats("tsv,vcf").unwrap();
        assert_eq!(fmts.len(), 2);
        assert_eq!(format_extension(fmts[0]), "tsv");
        assert_eq!(format_extension(fmts[1]), "vcf");
    }

    #[test]
    fn parse_output_formats_empty_defaults_to_tsv() {
        let fmts = parse_output_formats("").unwrap();
        assert_eq!(fmts.len(), 1);
        assert_eq!(format_extension(fmts[0]), "tsv");
    }

    #[test]
    fn parse_output_formats_unknown_errors() {
        let result = parse_output_formats("bam");
        assert!(result.is_err());
    }

    #[test]
    fn format_extension_all_variants() {
        assert_eq!(format_extension(OutputFormat::Tsv), "tsv");
        assert_eq!(format_extension(OutputFormat::Csv), "csv");
        assert_eq!(format_extension(OutputFormat::Json), "json");
        assert_eq!(format_extension(OutputFormat::Vcf), "vcf");
    }
}
