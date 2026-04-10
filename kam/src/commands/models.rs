//! `kam models` subcommand — list built-in ML models.

use crate::cli::{ModelsArgs, ModelsCommands};

/// Entry point for `kam models`.
///
/// # Errors
///
/// Returns an error if a subcommand fails.
pub fn run_models(args: ModelsArgs) -> Result<(), Box<dyn std::error::Error>> {
    match args.command {
        ModelsCommands::List => {
            println!("Built-in models:");
            for name in crate::models::builtin_names() {
                println!("  {name}");
            }
        }
    }
    Ok(())
}
