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
mod models;
mod output;
pub mod rescue;

fn main() {
    let cli = Cli::parse();

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
