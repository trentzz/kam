//! CLI definition for kam using clap derive API.
//!
//! # Example
//!
//! ```no_run
//! use clap::Parser;
//! use kam::cli::Cli;
//!
//! let cli = Cli::parse();
//! ```

use std::path::PathBuf;

use clap::{Parser, Subcommand};

use crate::caller_config::CallerConfigArgs;

/// Alignment-free variant detection for duplex UMI sequencing.
#[derive(Parser, Debug)]
#[command(
    name = "kam",
    about = "Alignment-free variant detection for duplex UMI sequencing"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

/// Available pipeline subcommands.
#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Assemble molecules from raw paired FASTQ input.
    Assemble(AssembleArgs),
    /// Build a k-mer index from consensus molecules against target sequences.
    Index(IndexArgs),
    /// Walk de Bruijn graph paths through indexed k-mers.
    Pathfind(PathfindArgs),
    /// Call variants from scored graph paths.
    Call(CallArgs),
    /// Run the full pipeline end-to-end.
    Run(Box<RunArgs>),
    /// List and manage built-in ML models.
    Models(ModelsArgs),
}

/// Arguments for the `assemble` subcommand.
#[derive(Parser, Debug)]
pub struct AssembleArgs {
    /// R1 input FASTQ file.
    #[arg(long)]
    pub r1: PathBuf,

    /// R2 input FASTQ file.
    #[arg(long)]
    pub r2: PathBuf,

    /// Output file for assembled molecules.
    #[arg(long)]
    pub output: PathBuf,

    /// Chemistry preset.
    #[arg(long, default_value = "twist-umi-duplex")]
    pub chemistry: String,

    /// Minimum Phred quality threshold for UMI bases.
    #[arg(long, default_value_t = 20u8)]
    pub min_umi_quality: u8,

    /// Minimum number of reads per UMI family.
    #[arg(long, default_value_t = 1u32)]
    pub min_family_size: u32,

    /// Minimum template length (optional).
    #[arg(long)]
    pub min_template_length: Option<u32>,

    /// Directory for log output.
    #[arg(long)]
    pub log_dir: Option<PathBuf>,

    /// Enable specific logs (repeatable).
    #[arg(long, action = clap::ArgAction::Append)]
    pub log: Vec<String>,

    /// Number of threads.
    #[arg(long)]
    pub threads: Option<usize>,
}

/// Arguments for the `index` subcommand.
#[derive(Parser, Debug)]
pub struct IndexArgs {
    /// Consensus molecules file (from `assemble`).
    #[arg(long)]
    pub input: PathBuf,

    /// Target sequences FASTA file.
    #[arg(long)]
    pub targets: PathBuf,

    /// Output k-mer index file.
    #[arg(long)]
    pub output: PathBuf,

    /// K-mer size.
    #[arg(short = 'k', long, default_value_t = 31u32)]
    pub kmer_size: u32,

    /// Optional FASTA of SV junction sequences to augment the k-mer allowlist.
    ///
    /// When provided, k-mers from these junction sequences are added to the
    /// allowlist alongside the target reference k-mers. This is required for
    /// detecting structural variants (deletions, tandem duplications, inversions)
    /// whose breakpoint k-mers are not in the reference target sequences.
    #[arg(long)]
    pub sv_junctions: Option<PathBuf>,
}

/// Arguments for the `pathfind` subcommand.
#[derive(Parser, Debug)]
pub struct PathfindArgs {
    /// K-mer index file (from `index`).
    #[arg(long)]
    pub index: PathBuf,

    /// Target sequences FASTA file.
    #[arg(long)]
    pub targets: PathBuf,

    /// Output scored paths file.
    #[arg(long)]
    pub output: PathBuf,

    /// K-mer size. Overrides the value inferred from the index when provided.
    ///
    /// By default, k is inferred from the stored k-mer values in the index.
    /// Use this flag when inference fails (e.g., all k-mers begin with A) or
    /// to explicitly confirm the expected k.
    #[arg(short = 'k', long)]
    pub kmer_size: Option<usize>,
}

/// Arguments for the `call` subcommand.
#[derive(Parser, Debug)]
pub struct CallArgs {
    /// Scored paths file (from `pathfind`).
    #[arg(long)]
    pub paths: PathBuf,

    /// Output variant calls file.
    #[arg(long)]
    pub output: PathBuf,

    /// Output format(s), comma-separated (tsv, csv, json, vcf).
    #[arg(long, default_value = "tsv")]
    pub output_format: String,

    /// Minimum posterior probability for a variant call.
    #[arg(long)]
    pub min_confidence: Option<f64>,

    /// Strand bias Fisher's exact p-value cutoff.
    #[arg(long)]
    pub strand_bias_threshold: Option<f64>,

    /// Minimum alt-supporting molecules to emit a call.
    #[arg(long)]
    pub min_alt_molecules: Option<u32>,

    /// Minimum variant-specific duplex molecules required for a PASS call.
    ///
    /// Counts only duplex molecules whose k-mers overlap the variant site
    /// (alt-path-specific k-mers).  Setting 1 requires at least one duplex
    /// confirmation on the alt allele.  Calls without sufficient duplex support
    /// are labelled LowDuplex.  Default: 0 (disabled).
    #[arg(long)]
    pub min_alt_duplex: Option<u32>,

    /// Minimum posterior probability for a structural variant PASS call.
    ///
    /// Applies to LargeDeletion, TandemDuplication, and Inversion types.
    /// Default: 0.95.
    #[arg(long)]
    pub sv_min_confidence: Option<f64>,

    /// Minimum alt-supporting molecules for a structural variant PASS call.
    ///
    /// Applies to LargeDeletion, TandemDuplication, and Inversion types.
    /// Default: 1.
    #[arg(long)]
    pub sv_min_alt_molecules: Option<u32>,

    /// Maximum VAF for a PASS call. Calls above this are labelled HighVaf.
    /// Useful for ctDNA somatic calling where germline heterozygous variants
    /// (VAF ≈ 0.5) should be excluded. Example: --max-vaf 0.35.
    #[arg(long)]
    pub max_vaf: Option<f64>,

    /// VCF file of expected somatic variants for tumour-informed monitoring mode.
    ///
    /// When provided, only calls whose (CHROM, POS, REF, ALT) matches an entry
    /// in this VCF are marked PASS. All other calls are labelled NotTargeted.
    /// Use this with a matched-tumour or reference-standard truth VCF to
    /// suppress background biological false positives.
    #[arg(long)]
    pub target_variants: Option<PathBuf>,

    /// Position tolerance (bp) for tumour-informed matching.
    ///
    /// When > 0, a call also passes the tumour-informed filter if its position
    /// is within this many bp of any target variant position, regardless of
    /// REF/ALT. Useful for large SVs where kam reports partial alleles that
    /// cannot match full SV truth alleles by exact REF/ALT. Default: 0
    /// (exact REF/ALT matching only).
    #[arg(long, default_value_t = 0u32)]
    pub ti_position_tolerance: u32,

    /// Fisher p-value threshold for strand bias filter on SV-type variants.
    ///
    /// Defaults to 1.0 (disabled). Inversion junction reads are structurally
    /// strand-biased and the standard threshold is inappropriate for SV paths.
    /// Set to 0.0 to apply the same threshold as SNVs/indels.
    #[arg(long, default_value_t = 1.0f64)]
    pub sv_strand_bias_threshold: f64,

    /// Built-in model name for optional variant re-scoring.
    ///
    /// Use a built-in name (e.g. `single-strand-v1`). The model is bundled
    /// with the binary and works out of the box after `cargo install`.
    /// The result is appended as `ml_prob` and `ml_filter` columns.
    /// Mutually exclusive with --custom-ml-model.
    #[arg(long)]
    pub ml_model: Option<String>,

    /// Path to a custom ONNX model file. The companion metadata file must
    /// exist at the same path with a `.json` extension.
    /// Mutually exclusive with --ml-model.
    #[arg(long)]
    pub custom_ml_model: Option<std::path::PathBuf>,
}

impl CallerConfigArgs for CallArgs {
    fn min_confidence(&self) -> Option<f64> {
        self.min_confidence
    }
    fn strand_bias_threshold(&self) -> Option<f64> {
        self.strand_bias_threshold
    }
    fn min_alt_molecules(&self) -> Option<u32> {
        self.min_alt_molecules
    }
    fn min_alt_duplex(&self) -> Option<u32> {
        self.min_alt_duplex
    }
    fn sv_min_confidence(&self) -> Option<f64> {
        self.sv_min_confidence
    }
    fn sv_min_alt_molecules(&self) -> Option<u32> {
        self.sv_min_alt_molecules
    }
    fn sv_strand_bias_threshold(&self) -> f64 {
        self.sv_strand_bias_threshold
    }
    fn max_vaf(&self) -> Option<f64> {
        self.max_vaf
    }
}

impl CallerConfigArgs for RunArgs {
    fn min_confidence(&self) -> Option<f64> {
        self.min_confidence
    }
    fn strand_bias_threshold(&self) -> Option<f64> {
        self.strand_bias_threshold
    }
    fn min_alt_molecules(&self) -> Option<u32> {
        self.min_alt_molecules
    }
    fn min_alt_duplex(&self) -> Option<u32> {
        self.min_alt_duplex
    }
    fn sv_min_confidence(&self) -> Option<f64> {
        self.sv_min_confidence
    }
    fn sv_min_alt_molecules(&self) -> Option<u32> {
        self.sv_min_alt_molecules
    }
    fn sv_strand_bias_threshold(&self) -> f64 {
        // Default to 1.0 (disabled) when not explicitly set.
        self.sv_strand_bias_threshold_override.unwrap_or(1.0)
    }
    fn max_vaf(&self) -> Option<f64> {
        self.max_vaf
    }
}

/// Arguments for the `run` subcommand (full pipeline).
///
/// All path inputs and most defaulted fields are `Option` so that a config
/// file (`--config`) can supply them.  When `--config` is absent, `r1`, `r2`,
/// `targets`, and `output_dir` are validated as required at runtime.
///
/// Priority: CLI flag > config file > built-in default.
#[derive(Parser, Debug)]
pub struct RunArgs {
    /// Path to a TOML config file.
    ///
    /// When provided, pipeline parameters are loaded from the file first, then
    /// any CLI flags override individual values.  When absent, all required
    /// fields must be supplied directly on the command line.
    #[arg(long)]
    pub config: Option<PathBuf>,

    /// R1 input FASTQ file.
    #[arg(long)]
    pub r1: Option<PathBuf>,

    /// R2 input FASTQ file.
    #[arg(long)]
    pub r2: Option<PathBuf>,

    /// Target sequences FASTA file.
    #[arg(long)]
    pub targets: Option<PathBuf>,

    /// Directory for all pipeline outputs.
    #[arg(long)]
    pub output_dir: Option<PathBuf>,

    /// Chemistry preset.
    ///
    /// Overrides the config file value when set.
    #[arg(long)]
    pub chemistry_override: Option<String>,

    /// Minimum Phred quality threshold for UMI bases.
    ///
    /// Overrides the config file value when set.
    #[arg(long)]
    pub min_umi_quality_override: Option<u8>,

    /// Minimum number of reads per UMI family.
    ///
    /// Overrides the config file value when set.
    #[arg(long)]
    pub min_family_size_override: Option<u32>,

    /// Minimum template length (optional).
    #[arg(long)]
    pub min_template_length: Option<u32>,

    /// K-mer size for indexing.
    ///
    /// Overrides the config file value when set.
    #[arg(short = 'k', long)]
    pub kmer_size_override: Option<u32>,

    /// Minimum posterior probability for a variant call.
    #[arg(long)]
    pub min_confidence: Option<f64>,

    /// Strand bias Fisher's exact p-value cutoff.
    #[arg(long)]
    pub strand_bias_threshold: Option<f64>,

    /// Minimum alt-supporting molecules to emit a call.
    #[arg(long)]
    pub min_alt_molecules: Option<u32>,

    /// Minimum posterior probability for a structural variant PASS call.
    ///
    /// Applies to LargeDeletion, TandemDuplication, and Inversion types.
    /// Default: 0.95.
    #[arg(long)]
    pub sv_min_confidence: Option<f64>,

    /// Minimum alt-supporting molecules for a structural variant PASS call.
    ///
    /// Applies to LargeDeletion, TandemDuplication, and Inversion types.
    /// Default: 1.
    #[arg(long)]
    pub sv_min_alt_molecules: Option<u32>,

    /// Minimum variant-specific duplex molecules required for a PASS call.
    ///
    /// Counts only duplex molecules whose k-mers overlap the variant site
    /// (alt-path-specific k-mers).  Setting 1 requires at least one duplex
    /// confirmation on the alt allele.  Calls without sufficient duplex support
    /// are labelled LowDuplex.  Default: 0 (disabled).
    #[arg(long)]
    pub min_alt_duplex: Option<u32>,

    /// Maximum VAF for a PASS call. Calls above this are labelled HighVaf.
    /// Useful for ctDNA somatic calling where germline heterozygous variants
    /// (VAF ≈ 0.5) should be excluded. Example: --max-vaf 0.35.
    #[arg(long)]
    pub max_vaf: Option<f64>,

    /// VCF file of expected somatic variants for tumour-informed monitoring mode.
    ///
    /// When provided, only calls whose (CHROM, POS, REF, ALT) matches an entry
    /// in this VCF are marked PASS. All other calls are labelled NotTargeted.
    /// Use this with a matched-tumour or reference-standard truth VCF to
    /// suppress background biological false positives.
    #[arg(long)]
    pub target_variants: Option<PathBuf>,

    /// Position tolerance (bp) for tumour-informed matching.
    ///
    /// When > 0, a call also passes the tumour-informed filter if its position
    /// is within this many bp of any target variant position, regardless of
    /// REF/ALT. Useful for large SVs where kam reports partial alleles that
    /// cannot match full SV truth alleles by exact REF/ALT. Default: 0
    /// (exact REF/ALT matching only).
    #[arg(long)]
    pub ti_position_tolerance_override: Option<u32>,

    /// Enable rescue probing for TI targets that produce no matching call.
    ///
    /// When set alongside `--target-variants`, the k-mer index is queried
    /// directly for each TI target variant that produces no PASS or
    /// sub-threshold call. Results appear in the output TSV with
    /// `call_source=RESCUED` or `call_source=NO_EVIDENCE`. Sub-threshold
    /// calls that match a TI target are marked `call_source=SUBTHRESHOLD`.
    #[arg(long, default_value_t = false)]
    pub ti_rescue: bool,

    /// Optional FASTA of SV junction sequences to augment the k-mer allowlist.
    ///
    /// When provided, k-mers from these junction sequences are added to the
    /// allowlist alongside the target reference k-mers. This is required for
    /// detecting structural variants (deletions, tandem duplications, inversions)
    /// whose breakpoint k-mers are not in the reference target sequences.
    #[arg(long)]
    pub sv_junctions: Option<PathBuf>,

    /// FASTA of synthetic fusion target sequences for fusion/translocation detection.
    ///
    /// Each entry must follow the format:
    /// `{name}__{chromA}:{startA}-{endA}__{chromB}:{startB}-{endB}__fusion`
    ///
    /// K-mers from these sequences are added to the allowlist so that
    /// fusion-spanning reads are captured. After normal target processing,
    /// fusion targets are walked and called separately using partner depth as
    /// the VAF denominator.
    #[arg(long)]
    pub fusion_targets: Option<PathBuf>,

    /// Output format(s), comma-separated (tsv, csv, json, vcf).
    ///
    /// Overrides the config file value when set.
    #[arg(long)]
    pub output_format_override: Option<String>,

    /// QC JSON output file.
    #[arg(long)]
    pub qc_output: Option<PathBuf>,

    /// Directory for log output.
    #[arg(long)]
    pub log_dir: Option<PathBuf>,

    /// Enable specific logs (repeatable).
    #[arg(long, action = clap::ArgAction::Append)]
    pub log: Vec<String>,

    /// Number of threads.
    #[arg(long)]
    pub threads: Option<usize>,

    /// Fisher p-value threshold for strand bias filter on SV-type variants.
    ///
    /// Overrides the config file value when set. Defaults to 1.0 (disabled).
    /// Inversion junction reads are structurally strand-biased and the standard
    /// threshold is inappropriate for SV paths.
    /// Set to 0.0 to apply the same threshold as SNVs/indels.
    #[arg(long)]
    pub sv_strand_bias_threshold_override: Option<f64>,

    /// Built-in model name for optional variant re-scoring.
    ///
    /// Use a built-in name (e.g. `single-strand-v1`). The model is bundled
    /// with the binary and works out of the box after `cargo install`.
    /// The result is appended as `ml_prob` and `ml_filter` columns.
    /// Mutually exclusive with --custom-ml-model.
    #[arg(long)]
    pub ml_model: Option<String>,

    /// Path to a custom ONNX model file. The companion metadata file must
    /// exist at the same path with a `.json` extension.
    /// Mutually exclusive with --ml-model.
    #[arg(long)]
    pub custom_ml_model: Option<std::path::PathBuf>,
}

/// Arguments for the `models` subcommand.
#[derive(clap::Args, Debug)]
pub struct ModelsArgs {
    #[command(subcommand)]
    pub command: ModelsCommands,
}

/// Subcommands under `kam models`.
#[derive(clap::Subcommand, Debug)]
pub enum ModelsCommands {
    /// List all built-in ML models.
    List,
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    /// Verify `--help` output can be generated without panicking.
    #[test]
    fn help_does_not_panic() {
        let mut cmd = Cli::command();
        let mut buf = Vec::new();
        cmd.write_help(&mut buf).expect("help should not fail");
        let help = String::from_utf8(buf).expect("help output should be valid UTF-8");
        assert!(help.contains("assemble"));
        assert!(help.contains("index"));
        assert!(help.contains("pathfind"));
        assert!(help.contains("call"));
        assert!(help.contains("run"));
    }

    /// Verify the `assemble` subcommand parses required arguments.
    #[test]
    fn assemble_parses_required_args() {
        let cli = Cli::try_parse_from([
            "kam",
            "assemble",
            "--r1",
            "r1.fq",
            "--r2",
            "r2.fq",
            "--output",
            "molecules.bin",
        ])
        .expect("assemble with required args should parse");
        let Commands::Assemble(args) = cli.command else {
            panic!("expected Assemble subcommand");
        };
        assert_eq!(args.r1, PathBuf::from("r1.fq"));
        assert_eq!(args.r2, PathBuf::from("r2.fq"));
        assert_eq!(args.output, PathBuf::from("molecules.bin"));
    }

    /// Verify `assemble` default values.
    #[test]
    fn assemble_defaults() {
        let cli = Cli::try_parse_from([
            "kam",
            "assemble",
            "--r1",
            "r1.fq",
            "--r2",
            "r2.fq",
            "--output",
            "molecules.bin",
        ])
        .expect("assemble defaults should parse");
        let Commands::Assemble(args) = cli.command else {
            panic!("expected Assemble subcommand");
        };
        assert_eq!(args.chemistry, "twist-umi-duplex");
        assert_eq!(args.min_umi_quality, 20);
        assert_eq!(args.min_family_size, 1);
        assert!(args.min_template_length.is_none());
        assert!(args.log_dir.is_none());
        assert!(args.threads.is_none());
        assert!(args.log.is_empty());
    }

    /// Verify the `index` subcommand parses required arguments.
    #[test]
    fn index_parses_required_args() {
        let cli = Cli::try_parse_from([
            "kam",
            "index",
            "--input",
            "molecules.bin",
            "--targets",
            "targets.fa",
            "--output",
            "index.bin",
        ])
        .expect("index with required args should parse");
        let Commands::Index(args) = cli.command else {
            panic!("expected Index subcommand");
        };
        assert_eq!(args.input, PathBuf::from("molecules.bin"));
        assert_eq!(args.targets, PathBuf::from("targets.fa"));
        assert_eq!(args.output, PathBuf::from("index.bin"));
    }

    /// Verify `index` default k-mer size.
    #[test]
    fn index_defaults() {
        let cli = Cli::try_parse_from([
            "kam",
            "index",
            "--input",
            "molecules.bin",
            "--targets",
            "targets.fa",
            "--output",
            "index.bin",
        ])
        .expect("index defaults should parse");
        let Commands::Index(args) = cli.command else {
            panic!("expected Index subcommand");
        };
        assert_eq!(args.kmer_size, 31);
    }

    /// Verify the `pathfind` subcommand parses required arguments.
    #[test]
    fn pathfind_parses_required_args() {
        let cli = Cli::try_parse_from([
            "kam",
            "pathfind",
            "--index",
            "index.bin",
            "--targets",
            "targets.fa",
            "--output",
            "paths.bin",
        ])
        .expect("pathfind with required args should parse");
        let Commands::Pathfind(args) = cli.command else {
            panic!("expected Pathfind subcommand");
        };
        assert_eq!(args.index, PathBuf::from("index.bin"));
        assert_eq!(args.targets, PathBuf::from("targets.fa"));
        assert_eq!(args.output, PathBuf::from("paths.bin"));
    }

    /// Verify the `call` subcommand parses required arguments.
    #[test]
    fn call_parses_required_args() {
        let cli = Cli::try_parse_from([
            "kam",
            "call",
            "--paths",
            "paths.bin",
            "--output",
            "variants.tsv",
        ])
        .expect("call with required args should parse");
        let Commands::Call(args) = cli.command else {
            panic!("expected Call subcommand");
        };
        assert_eq!(args.paths, PathBuf::from("paths.bin"));
        assert_eq!(args.output, PathBuf::from("variants.tsv"));
    }

    /// Verify `call` default output format.
    #[test]
    fn call_defaults() {
        let cli = Cli::try_parse_from([
            "kam",
            "call",
            "--paths",
            "paths.bin",
            "--output",
            "variants.tsv",
        ])
        .expect("call defaults should parse");
        let Commands::Call(args) = cli.command else {
            panic!("expected Call subcommand");
        };
        assert_eq!(args.output_format, "tsv");
        assert!(args.min_confidence.is_none());
        assert!(args.strand_bias_threshold.is_none());
        assert_eq!(args.sv_strand_bias_threshold, 1.0);
    }

    /// Verify the `run` subcommand parses CLI path arguments.
    #[test]
    fn run_parses_required_args() {
        let cli = Cli::try_parse_from([
            "kam",
            "run",
            "--r1",
            "r1.fq",
            "--r2",
            "r2.fq",
            "--targets",
            "targets.fa",
            "--output-dir",
            "outdir",
        ])
        .expect("run with required args should parse");
        let Commands::Run(args) = cli.command else {
            panic!("expected Run subcommand");
        };
        assert_eq!(args.r1, Some(PathBuf::from("r1.fq")));
        assert_eq!(args.r2, Some(PathBuf::from("r2.fq")));
        assert_eq!(args.targets, Some(PathBuf::from("targets.fa")));
        assert_eq!(args.output_dir, Some(PathBuf::from("outdir")));
    }

    /// Verify `run` optional override fields are None by default.
    #[test]
    fn run_defaults() {
        let cli = Cli::try_parse_from([
            "kam",
            "run",
            "--r1",
            "r1.fq",
            "--r2",
            "r2.fq",
            "--targets",
            "targets.fa",
            "--output-dir",
            "outdir",
        ])
        .expect("run defaults should parse");
        let Commands::Run(args) = cli.command else {
            panic!("expected Run subcommand");
        };
        // Override fields are None when not explicitly provided.
        assert!(args.chemistry_override.is_none());
        assert!(args.min_umi_quality_override.is_none());
        assert!(args.min_family_size_override.is_none());
        assert!(args.kmer_size_override.is_none());
        assert!(args.output_format_override.is_none());
        assert!(args.qc_output.is_none());
        assert!(args.sv_strand_bias_threshold_override.is_none());
    }

    /// Verify `--config` flag is accepted.
    #[test]
    fn run_config_flag_accepted() {
        let cli = Cli::try_parse_from(["kam", "run", "--config", "config.toml"])
            .expect("run --config should parse");
        let Commands::Run(args) = cli.command else {
            panic!("expected Run subcommand");
        };
        assert_eq!(args.config, Some(PathBuf::from("config.toml")));
    }

    /// Verify override flags are captured when supplied.
    #[test]
    fn run_override_flags_captured() {
        let cli = Cli::try_parse_from([
            "kam",
            "run",
            "--r1",
            "r1.fq",
            "--chemistry-override",
            "custom",
            "--kmer-size-override",
            "21",
            "--output-format-override",
            "vcf",
        ])
        .expect("run override flags should parse");
        let Commands::Run(args) = cli.command else {
            panic!("expected Run subcommand");
        };
        assert_eq!(args.chemistry_override, Some("custom".to_string()));
        assert_eq!(args.kmer_size_override, Some(21u32));
        assert_eq!(args.output_format_override, Some("vcf".to_string()));
    }

    /// Verify missing required arguments produce an error (not a panic).
    #[test]
    fn assemble_missing_required_arg_errors() {
        let result = Cli::try_parse_from(["kam", "assemble", "--r1", "r1.fq"]);
        assert!(
            result.is_err(),
            "missing --r2 and --output should be an error"
        );
    }

    /// Verify repeatable --log flag accumulates values.
    #[test]
    fn assemble_log_flag_repeatable() {
        let cli = Cli::try_parse_from([
            "kam",
            "assemble",
            "--r1",
            "r1.fq",
            "--r2",
            "r2.fq",
            "--output",
            "molecules.bin",
            "--log",
            "umi",
            "--log",
            "family",
        ])
        .expect("repeatable --log should parse");
        let Commands::Assemble(args) = cli.command else {
            panic!("expected Assemble subcommand");
        };
        assert_eq!(args.log, vec!["umi", "family"]);
    }
}
