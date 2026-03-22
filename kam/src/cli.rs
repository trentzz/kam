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
}

/// Arguments for the `run` subcommand (full pipeline).
#[derive(Parser, Debug)]
pub struct RunArgs {
    /// R1 input FASTQ file.
    #[arg(long)]
    pub r1: PathBuf,

    /// R2 input FASTQ file.
    #[arg(long)]
    pub r2: PathBuf,

    /// Target sequences FASTA file.
    #[arg(long)]
    pub targets: PathBuf,

    /// Directory for all pipeline outputs.
    #[arg(long)]
    pub output_dir: PathBuf,

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

    /// K-mer size for indexing.
    #[arg(short = 'k', long, default_value_t = 31u32)]
    pub kmer_size: u32,

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

    /// Optional FASTA of SV junction sequences to augment the k-mer allowlist.
    ///
    /// When provided, k-mers from these junction sequences are added to the
    /// allowlist alongside the target reference k-mers. This is required for
    /// detecting structural variants (deletions, tandem duplications, inversions)
    /// whose breakpoint k-mers are not in the reference target sequences.
    #[arg(long)]
    pub sv_junctions: Option<PathBuf>,

    /// Output format(s), comma-separated (tsv, csv, json, vcf).
    #[arg(long, default_value = "tsv")]
    pub output_format: String,

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
    }

    /// Verify the `run` subcommand parses required arguments.
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
        assert_eq!(args.r1, PathBuf::from("r1.fq"));
        assert_eq!(args.r2, PathBuf::from("r2.fq"));
        assert_eq!(args.targets, PathBuf::from("targets.fa"));
        assert_eq!(args.output_dir, PathBuf::from("outdir"));
    }

    /// Verify `run` default values.
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
        assert_eq!(args.chemistry, "twist-umi-duplex");
        assert_eq!(args.min_umi_quality, 20);
        assert_eq!(args.min_family_size, 1);
        assert_eq!(args.kmer_size, 31);
        assert_eq!(args.output_format, "tsv");
        assert!(args.qc_output.is_none());
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
