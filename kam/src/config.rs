//! Unified configuration file support for `kam run`.
//!
//! Provides [`KamConfig`], a TOML-serialisable struct that covers all
//! pipeline parameters.  CLI flags always take precedence over config file
//! values, which in turn take precedence over hard-coded defaults.
//!
//! # Example
//!
//! ```no_run
//! use kam::config::KamConfig;
//!
//! let cfg = KamConfig::from_file("config.toml").expect("load config");
//! ```

use std::path::PathBuf;

use serde::Deserialize;

use kam_assemble::assembler::AssemblerConfig;
use kam_assemble::consensus::ConsensusConfig;
use kam_assemble::parser::ParserConfig;
use kam_call::caller::CallerConfig;
use kam_core::chemistry::ReadStructure;

use crate::cli::RunArgs;

// ── Section structs ───────────────────────────────────────────────────────────

/// Input file paths.
///
/// All fields are optional here so a partial config file is still valid.
/// The [`KamConfig::validate`] method checks that required paths are present
/// before the pipeline starts.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InputConfig {
    /// R1 FASTQ input file.
    pub r1: Option<PathBuf>,
    /// R2 FASTQ input file.
    pub r2: Option<PathBuf>,
    /// Target sequences FASTA file.
    pub targets: Option<PathBuf>,
    /// Optional SV junction sequences FASTA.
    pub sv_junctions: Option<PathBuf>,
    /// VCF of expected somatic variants for tumour-informed monitoring mode.
    pub target_variants: Option<PathBuf>,
    /// FASTA of synthetic fusion target sequences.
    ///
    /// Each entry encodes both partner coordinates in its ID using the format:
    /// `{name}__{chromA}:{startA}-{endA}__{chromB}:{startB}-{endB}__fusion`.
    /// K-mers from these sequences are added to the allowlist alongside normal
    /// target k-mers so that fusion-spanning reads are captured during indexing.
    pub fusion_targets: Option<PathBuf>,
    /// FASTA of raw junction sequences for allowlist augmentation and standalone walking.
    ///
    /// Each sequence is added to the k-mer allowlist and walked as a standalone
    /// target using total library depth as the VAF denominator. Useful when a
    /// junction sequence is observed directly in a BAM or IGV session without
    /// known genomic coordinates.
    pub junction_sequences: Option<PathBuf>,
}

/// Output settings.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OutputConfig {
    /// Directory for all pipeline outputs.
    pub output_dir: Option<PathBuf>,
    /// Output format string, comma-separated (tsv, csv, json, vcf).
    pub output_format: Option<String>,
    /// Optional QC JSON output path.
    pub qc_output: Option<PathBuf>,
}

impl Default for OutputConfig {
    fn default() -> Self {
        Self {
            output_dir: None,
            output_format: Some("tsv".to_string()),
            qc_output: None,
        }
    }
}

/// Chemistry and read-structure settings.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ChemistryConfig {
    /// Chemistry preset name.
    pub preset: Option<String>,
    /// UMI length in bases (default: 5 for Twist).
    #[serde(default = "default_umi_length")]
    pub umi_length: u32,
    /// Skip (spacer) length in bases (default: 2 for Twist).
    #[serde(default = "default_skip_length")]
    pub skip_length: u32,
    /// Whether the chemistry is duplex (default: true).
    #[serde(default = "default_duplex")]
    pub duplex: bool,
    /// Minimum Phred quality for UMI bases.
    pub min_umi_quality: Option<u8>,
    /// Minimum template length in bases.
    pub min_template_length: Option<u32>,
}

fn default_umi_length() -> u32 {
    5
}
fn default_skip_length() -> u32 {
    2
}
fn default_duplex() -> bool {
    true
}

impl Default for ChemistryConfig {
    fn default() -> Self {
        Self {
            preset: Some("twist-umi-duplex".to_string()),
            umi_length: 5,
            skip_length: 2,
            duplex: true,
            min_umi_quality: Some(20),
            min_template_length: None,
        }
    }
}

/// Molecule assembly settings.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AssemblyConfig {
    /// Minimum reads per UMI family.
    pub min_family_size: Option<u32>,
}

impl Default for AssemblyConfig {
    fn default() -> Self {
        Self {
            min_family_size: Some(1),
        }
    }
}

/// K-mer indexing settings.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IndexingConfig {
    /// K-mer size.
    pub kmer_size: Option<u32>,
}

impl Default for IndexingConfig {
    fn default() -> Self {
        Self {
            kmer_size: Some(31),
        }
    }
}

/// Variant calling thresholds.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CallingConfig {
    /// Minimum posterior probability to emit a PASS call.
    pub min_confidence: Option<f64>,
    /// Fisher p-value cutoff for strand bias.
    pub strand_bias_threshold: Option<f64>,
    /// Minimum alt-supporting molecules required.
    pub min_alt_molecules: Option<u32>,
    /// Minimum variant-specific duplex molecules required.
    pub min_alt_duplex: Option<u32>,
    /// Maximum VAF for a PASS call.
    pub max_vaf: Option<f64>,
    /// Minimum posterior probability for SV-type calls.
    pub sv_min_confidence: Option<f64>,
    /// Minimum alt molecules for SV-type calls.
    pub sv_min_alt_molecules: Option<u32>,
    /// Strand bias threshold for SV-type variants.
    pub sv_strand_bias_threshold: Option<f64>,
    /// Position tolerance (bp) for tumour-informed matching.
    pub ti_position_tolerance: Option<u32>,
    /// Enable k-mer rescue probe for undetected TI targets.
    #[serde(default)]
    pub ti_rescue: bool,
}

/// Logging settings.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LoggingConfig {
    /// Directory for log output.
    pub log_dir: Option<PathBuf>,
    /// Specific log channels to enable.
    pub log: Option<Vec<String>>,
}

/// Runtime settings.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RuntimeConfig {
    /// Number of threads.
    pub threads: Option<usize>,
}

// ── Top-level config ──────────────────────────────────────────────────────────

/// Unified configuration for the `kam run` pipeline.
///
/// Deserialised from a TOML file, then merged with any CLI overrides.
/// CLI flags always win over config file values.
///
/// # Example
///
/// ```
/// use kam::config::KamConfig;
///
/// let toml = r#"
/// [input]
/// r1 = "R1.fq"
/// r2 = "R2.fq"
/// targets = "targets.fa"
///
/// [output]
/// output_dir = "results"
/// "#;
///
/// let cfg: KamConfig = toml::from_str(toml).expect("parse");
/// assert_eq!(cfg.input.r1, Some("R1.fq".into()));
/// ```
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct KamConfig {
    /// Input file paths.
    #[serde(default)]
    pub input: InputConfig,
    /// Output settings.
    #[serde(default)]
    pub output: OutputConfig,
    /// Chemistry and read-structure settings.
    #[serde(default)]
    pub chemistry: ChemistryConfig,
    /// Molecule assembly settings.
    #[serde(default)]
    pub assembly: AssemblyConfig,
    /// K-mer indexing settings.
    #[serde(default)]
    pub indexing: IndexingConfig,
    /// Variant calling thresholds.
    #[serde(default)]
    pub calling: CallingConfig,
    /// Logging settings.
    #[serde(default)]
    pub logging: LoggingConfig,
    /// Runtime settings.
    #[serde(default)]
    pub runtime: RuntimeConfig,
}

impl KamConfig {
    /// Load a [`KamConfig`] from a TOML file at `path`.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be read or the TOML is malformed.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use kam::config::KamConfig;
    ///
    /// let cfg = KamConfig::from_file("config.toml").expect("load config");
    /// ```
    pub fn from_file(
        path: impl AsRef<std::path::Path>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let text = std::fs::read_to_string(path)?;
        let cfg: Self = toml::from_str(&text)?;
        Ok(cfg)
    }

    /// Build a [`KamConfig`] directly from [`RunArgs`].
    ///
    /// Used when no `--config` file is given. Every CLI flag maps to the
    /// corresponding config field; the result is functionally identical to the
    /// original CLI-only behaviour.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use std::path::PathBuf;
    /// use kam::config::KamConfig;
    /// use kam::cli::RunArgs;
    ///
    /// // Build via clap in real usage; constructed directly here for demonstration.
    /// ```
    pub fn from_cli(args: &RunArgs) -> Self {
        let mut cfg = Self::default();

        cfg.input.r1 = args.r1.clone();
        cfg.input.r2 = args.r2.clone();
        cfg.input.targets = args.targets.clone();
        cfg.input.sv_junctions = args.sv_junctions.clone();
        cfg.input.target_variants = args.target_variants.clone();
        cfg.input.fusion_targets = args.fusion_targets.clone();
        cfg.input.junction_sequences = args.junction_sequences.clone();

        cfg.output.output_dir = args.output_dir.clone();
        cfg.output.output_format = args.output_format_override.clone();
        cfg.output.qc_output = args.qc_output.clone();

        cfg.chemistry.preset = args.chemistry_override.clone();
        cfg.chemistry.min_umi_quality = args.min_umi_quality_override;
        cfg.chemistry.min_template_length = args.min_template_length;

        cfg.assembly.min_family_size = args.min_family_size_override;

        cfg.indexing.kmer_size = args.kmer_size_override;

        cfg.calling.min_confidence = args.min_confidence;
        cfg.calling.strand_bias_threshold = args.strand_bias_threshold;
        cfg.calling.min_alt_molecules = args.min_alt_molecules;
        cfg.calling.min_alt_duplex = args.min_alt_duplex;
        cfg.calling.max_vaf = args.max_vaf;
        cfg.calling.sv_min_confidence = args.sv_min_confidence;
        cfg.calling.sv_min_alt_molecules = args.sv_min_alt_molecules;
        cfg.calling.sv_strand_bias_threshold = args.sv_strand_bias_threshold_override;
        cfg.calling.ti_position_tolerance = args.ti_position_tolerance_override;
        cfg.calling.ti_rescue = args.ti_rescue;

        cfg.logging.log_dir = args.log_dir.clone();
        cfg.logging.log = if args.log.is_empty() {
            None
        } else {
            Some(args.log.clone())
        };

        cfg.runtime.threads = args.threads;

        cfg
    }

    /// Apply CLI overrides on top of an existing config.
    ///
    /// Any CLI flag that was explicitly provided (i.e. is `Some`) wins over
    /// the current config value. For non-`Option` fields with defaults (e.g.
    /// `--chemistry`, `--kmer-size`), we compare against the clap default and
    /// only override when the user passed the flag explicitly — but in
    /// practice the simplest safe rule is to always apply CLI non-option
    /// fields, because the user has to supply them on the command line anyway.
    ///
    /// Required fields (`r1`, `r2`, `targets`, `output_dir`) use the same
    /// pattern: `Some` from CLI overrides whatever the config file says.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use kam::config::KamConfig;
    ///
    /// let mut cfg = KamConfig::default();
    /// // args would come from clap in real usage
    /// ```
    pub fn merge_cli_overrides(&mut self, args: &RunArgs) {
        // Required path fields: only override when CLI value is Some.
        if let Some(ref v) = args.r1 {
            self.input.r1 = Some(v.clone());
        }
        if let Some(ref v) = args.r2 {
            self.input.r2 = Some(v.clone());
        }
        if let Some(ref v) = args.targets {
            self.input.targets = Some(v.clone());
        }
        if let Some(ref v) = args.sv_junctions {
            self.input.sv_junctions = Some(v.clone());
        }
        if let Some(ref v) = args.target_variants {
            self.input.target_variants = Some(v.clone());
        }
        if let Some(ref v) = args.fusion_targets {
            self.input.fusion_targets = Some(v.clone());
        }
        if let Some(ref v) = args.junction_sequences {
            self.input.junction_sequences = Some(v.clone());
        }

        // Output fields.
        if let Some(ref v) = args.output_dir {
            self.output.output_dir = Some(v.clone());
        }
        if let Some(ref v) = &args.output_format_override {
            self.output.output_format = Some(v.clone());
        }
        if let Some(ref v) = args.qc_output {
            self.output.qc_output = Some(v.clone());
        }

        // Chemistry fields.
        if let Some(ref v) = args.chemistry_override {
            self.chemistry.preset = Some(v.clone());
        }
        if let Some(v) = args.min_umi_quality_override {
            self.chemistry.min_umi_quality = Some(v);
        }
        if let Some(v) = args.min_template_length {
            self.chemistry.min_template_length = Some(v);
        }

        // Assembly fields.
        if let Some(v) = args.min_family_size_override {
            self.assembly.min_family_size = Some(v);
        }

        // Indexing fields.
        if let Some(v) = args.kmer_size_override {
            self.indexing.kmer_size = Some(v);
        }

        // Calling fields.
        if let Some(v) = args.min_confidence {
            self.calling.min_confidence = Some(v);
        }
        if let Some(v) = args.strand_bias_threshold {
            self.calling.strand_bias_threshold = Some(v);
        }
        if let Some(v) = args.min_alt_molecules {
            self.calling.min_alt_molecules = Some(v);
        }
        if let Some(v) = args.min_alt_duplex {
            self.calling.min_alt_duplex = Some(v);
        }
        if let Some(v) = args.max_vaf {
            self.calling.max_vaf = Some(v);
        }
        if let Some(v) = args.sv_min_confidence {
            self.calling.sv_min_confidence = Some(v);
        }
        if let Some(v) = args.sv_min_alt_molecules {
            self.calling.sv_min_alt_molecules = Some(v);
        }
        if let Some(v) = args.sv_strand_bias_threshold_override {
            self.calling.sv_strand_bias_threshold = Some(v);
        }
        if let Some(v) = args.ti_position_tolerance_override {
            self.calling.ti_position_tolerance = Some(v);
        }
        if args.ti_rescue {
            self.calling.ti_rescue = true;
        }

        // Logging fields.
        if let Some(ref v) = args.log_dir {
            self.logging.log_dir = Some(v.clone());
        }
        if !args.log.is_empty() {
            self.logging.log = Some(args.log.clone());
        }

        // Runtime fields.
        if let Some(v) = args.threads {
            self.runtime.threads = Some(v);
        }
    }

    /// Validate that all required fields are present.
    ///
    /// Returns an error message naming the first missing required field.
    ///
    /// # Errors
    ///
    /// Returns `Err` if `r1`, `r2`, `targets`, or `output_dir` is absent.
    ///
    /// # Example
    ///
    /// ```
    /// use kam::config::KamConfig;
    ///
    /// let mut cfg = KamConfig::default();
    /// assert!(cfg.validate().is_err());
    /// ```
    pub fn validate(&self) -> Result<(), String> {
        if self.input.r1.is_none() {
            return Err(
                "missing required field: r1 (provide --r1 or set [input] r1 in config)".to_string(),
            );
        }
        if self.input.r2.is_none() {
            return Err(
                "missing required field: r2 (provide --r2 or set [input] r2 in config)".to_string(),
            );
        }
        if self.input.targets.is_none() {
            return Err("missing required field: targets (provide --targets or set [input] targets in config)".to_string());
        }
        if self.output.output_dir.is_none() {
            return Err("missing required field: output_dir (provide --output-dir or set [output] output_dir in config)".to_string());
        }
        Ok(())
    }

    /// Convert this config into a [`ParserConfig`].
    ///
    /// # Example
    ///
    /// ```
    /// use kam::config::KamConfig;
    ///
    /// let cfg = KamConfig::default();
    /// let parser = cfg.to_parser_config();
    /// assert!(parser.min_umi_quality.is_some());
    /// ```
    pub fn to_parser_config(&self) -> ParserConfig {
        ParserConfig {
            read_structure: ReadStructure {
                umi_length: self.chemistry.umi_length as usize,
                skip_length: self.chemistry.skip_length as usize,
            },
            min_template_length: self.chemistry.min_template_length.map(|v| v as usize),
            min_umi_quality: self.chemistry.min_umi_quality.and_then(|q| {
                if q == 0 {
                    None
                } else {
                    Some(q)
                }
            }),
        }
    }

    /// Convert this config into an [`AssemblerConfig`].
    ///
    /// # Example
    ///
    /// ```
    /// use kam::config::KamConfig;
    ///
    /// let cfg = KamConfig::default();
    /// let asm = cfg.to_assembler_config();
    /// assert_eq!(asm.min_family_size, 1);
    /// ```
    pub fn to_assembler_config(&self) -> AssemblerConfig {
        AssemblerConfig {
            min_family_size: self.assembly.min_family_size.unwrap_or(1) as u8,
            consensus: ConsensusConfig::default(),
            ..AssemblerConfig::default()
        }
    }

    /// Convert this config into a [`CallerConfig`].
    ///
    /// Falls back to [`CallerConfig::default`] for any field not set.
    ///
    /// # Example
    ///
    /// ```
    /// use kam::config::KamConfig;
    ///
    /// let cfg = KamConfig::default();
    /// let caller = cfg.to_caller_config();
    /// assert_eq!(caller.min_confidence, 0.99);
    /// ```
    pub fn to_caller_config(&self) -> CallerConfig {
        let defaults = CallerConfig::default();
        CallerConfig {
            min_confidence: self
                .calling
                .min_confidence
                .unwrap_or(defaults.min_confidence),
            strand_bias_threshold: self
                .calling
                .strand_bias_threshold
                .unwrap_or(defaults.strand_bias_threshold),
            min_alt_molecules: self
                .calling
                .min_alt_molecules
                .unwrap_or(defaults.min_alt_molecules),
            min_alt_duplex: self
                .calling
                .min_alt_duplex
                .unwrap_or(defaults.min_alt_duplex),
            max_vaf: self.calling.max_vaf.or(defaults.max_vaf),
            sv_min_confidence: self
                .calling
                .sv_min_confidence
                .unwrap_or(defaults.sv_min_confidence),
            sv_min_alt_molecules: self
                .calling
                .sv_min_alt_molecules
                .unwrap_or(defaults.sv_min_alt_molecules),
            sv_strand_bias_threshold: self
                .calling
                .sv_strand_bias_threshold
                .unwrap_or(defaults.sv_strand_bias_threshold),
            ..defaults
        }
    }

    /// Return the k-mer size, falling back to the default of 31.
    ///
    /// # Example
    ///
    /// ```
    /// use kam::config::KamConfig;
    ///
    /// let cfg = KamConfig::default();
    /// assert_eq!(cfg.kmer_size(), 31);
    /// ```
    pub fn kmer_size(&self) -> usize {
        self.indexing.kmer_size.unwrap_or(31) as usize
    }

    /// Return the output format string, falling back to `"tsv"`.
    ///
    /// # Example
    ///
    /// ```
    /// use kam::config::KamConfig;
    ///
    /// let cfg = KamConfig::default();
    /// assert_eq!(cfg.output_format(), "tsv");
    /// ```
    pub fn output_format(&self) -> &str {
        self.output.output_format.as_deref().unwrap_or("tsv")
    }

    /// Return the `ti_position_tolerance`, falling back to 0.
    ///
    /// # Example
    ///
    /// ```
    /// use kam::config::KamConfig;
    ///
    /// let cfg = KamConfig::default();
    /// assert_eq!(cfg.ti_position_tolerance(), 0);
    /// ```
    pub fn ti_position_tolerance(&self) -> i64 {
        self.calling.ti_position_tolerance.unwrap_or(0) as i64
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn minimal_config_str() -> &'static str {
        r#"
[input]
r1 = "R1.fq"
r2 = "R2.fq"
targets = "targets.fa"

[output]
output_dir = "results"
"#
    }

    /// Minimal TOML parses correctly.
    #[test]
    fn parse_minimal_toml() {
        let cfg: KamConfig = toml::from_str(minimal_config_str()).expect("parse");
        assert_eq!(cfg.input.r1, Some(PathBuf::from("R1.fq")));
        assert_eq!(cfg.input.r2, Some(PathBuf::from("R2.fq")));
        assert_eq!(cfg.input.targets, Some(PathBuf::from("targets.fa")));
        assert_eq!(cfg.output.output_dir, Some(PathBuf::from("results")));
    }

    /// Validate passes when all required fields are present.
    #[test]
    fn validate_passes_when_complete() {
        let cfg: KamConfig = toml::from_str(minimal_config_str()).expect("parse");
        assert!(cfg.validate().is_ok());
    }

    /// Validate fails when r1 is missing.
    #[test]
    fn validate_fails_missing_r1() {
        let toml = r#"
[input]
r2 = "R2.fq"
targets = "targets.fa"
[output]
output_dir = "results"
"#;
        let cfg: KamConfig = toml::from_str(toml).expect("parse");
        assert!(cfg.validate().is_err());
    }

    /// Validate fails when output_dir is missing.
    #[test]
    fn validate_fails_missing_output_dir() {
        let cfg: KamConfig = toml::from_str(minimal_config_str())
            .map(|mut c: KamConfig| {
                c.output.output_dir = None;
                c
            })
            .expect("parse");
        assert!(cfg.validate().is_err());
    }

    /// Config file values are used when CLI flags are absent.
    #[test]
    fn config_file_values_used_when_cli_absent() {
        let toml = r#"
[input]
r1 = "file_r1.fq"
r2 = "file_r2.fq"
targets = "file_targets.fa"
[output]
output_dir = "file_out"
[calling]
min_confidence = 0.95
"#;
        let cfg: KamConfig = toml::from_str(toml).expect("parse");
        assert_eq!(cfg.calling.min_confidence, Some(0.95));
        let caller = cfg.to_caller_config();
        assert_eq!(caller.min_confidence, 0.95);
    }

    /// to_parser_config maps chemistry fields correctly.
    #[test]
    fn to_parser_config_maps_fields() {
        let toml = r#"
[input]
r1 = "r1.fq"
r2 = "r2.fq"
targets = "t.fa"
[output]
output_dir = "out"
[chemistry]
min_umi_quality = 25
min_template_length = 50
"#;
        let cfg: KamConfig = toml::from_str(toml).expect("parse");
        let parser = cfg.to_parser_config();
        assert_eq!(parser.min_umi_quality, Some(25));
        assert_eq!(parser.min_template_length, Some(50));
    }

    /// to_assembler_config maps assembly fields correctly.
    #[test]
    fn to_assembler_config_maps_fields() {
        let toml = r#"
[input]
r1 = "r1.fq"
r2 = "r2.fq"
targets = "t.fa"
[output]
output_dir = "out"
[assembly]
min_family_size = 3
"#;
        let cfg: KamConfig = toml::from_str(toml).expect("parse");
        let asm = cfg.to_assembler_config();
        assert_eq!(asm.min_family_size, 3);
    }

    /// to_caller_config falls back to CallerConfig defaults for unset fields.
    #[test]
    fn to_caller_config_defaults_for_unset() {
        let cfg: KamConfig = toml::from_str(minimal_config_str()).expect("parse");
        let caller = cfg.to_caller_config();
        let defaults = CallerConfig::default();
        assert_eq!(caller.min_confidence, defaults.min_confidence);
        assert_eq!(caller.strand_bias_threshold, defaults.strand_bias_threshold);
        assert_eq!(caller.min_alt_molecules, defaults.min_alt_molecules);
    }

    /// kmer_size falls back to 31 when not set.
    #[test]
    fn kmer_size_defaults_to_31() {
        let cfg: KamConfig = toml::from_str(minimal_config_str()).expect("parse");
        assert_eq!(cfg.kmer_size(), 31);
    }

    /// output_format falls back to "tsv" when not set.
    #[test]
    fn output_format_defaults_to_tsv() {
        let cfg: KamConfig = toml::from_str(minimal_config_str()).expect("parse");
        assert_eq!(cfg.output_format(), "tsv");
    }
}
