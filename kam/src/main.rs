//! kam — alignment-free variant detection for duplex UMI sequencing.
//!
//! Integrates kam-core, kam-assemble, kam-index, kam-pathfind, and kam-call into a
//! single binary with a Nextflow-friendly CLI.

use clap::Parser;

use crate::cli::{Cli, Commands};

mod caller_config;
mod cli;
mod commands;
pub mod config;
pub mod memory_budget;
pub mod metrics;
mod models;
mod output;
pub mod rescue;

fn main() {
    let cli = Cli::parse();

    // Determine log level from CLI --log-level flag on subcommands that support it.
    let log_level = match &cli.command {
        Commands::Run(args) => args.log_level.as_deref(),
        Commands::Assemble(args) => args.log_level.as_deref(),
        _ => None,
    };

    // Initialise env_logger. CLI --log-level takes precedence, then RUST_LOG
    // env var, then defaults to "info".
    let mut builder =
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"));
    if let Some(level) = log_level {
        builder.parse_filters(level);
    }
    builder.init();

    let result = match cli.command {
        Commands::Assemble(args) => commands::assemble::run_assemble(args),
        Commands::Index(args) => commands::index::run_index(args),
        Commands::Pathfind(args) => commands::pathfind::run_pathfind(args),
        Commands::Call(args) => commands::call::run_call(args),
        Commands::Run(args) => commands::run::run_pipeline(*args),
        Commands::Models(args) => commands::models::run_models(args),
        Commands::Explore(args) => commands::explore::run_explore(args),
    };

    if let Err(e) = result {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}
