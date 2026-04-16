//! `kam explore` subcommand — interactive REPL for querying bincode data files.
//!
//! Loads a kam bincode file (molecules, k-mer index, scored paths), then
//! presents an interactive prompt where users can inspect, filter, and export
//! records without writing one-off scripts.

use std::fs;
use std::io::Write;
use std::path::Path;

use comfy_table::{Cell, Table};
use rustyline::error::ReadlineError;

use kam_core::molecule::{FamilyType, Molecule};
use kam_core::serialize::{read_bincode, read_header, FileType};

use crate::cli::ExploreArgs;
use crate::commands::index::KmerEntry;
use crate::commands::pathfind::ScoredPathRecord;

// ─── Data container ─────────────────────────────────────────────────────────

/// Union of supported bincode record types.
enum ExploreData {
    Molecules(Vec<Molecule>),
    KmerIndex(Vec<KmerEntry>),
    ScoredPaths(Vec<ScoredPathRecord>),
}

impl ExploreData {
    /// Number of records loaded.
    fn len(&self) -> usize {
        match self {
            ExploreData::Molecules(v) => v.len(),
            ExploreData::KmerIndex(v) => v.len(),
            ExploreData::ScoredPaths(v) => v.len(),
        }
    }

    /// Human-readable label for the file type.
    fn type_label(&self) -> &'static str {
        match self {
            ExploreData::Molecules(_) => "Molecules",
            ExploreData::KmerIndex(_) => "KmerIndex",
            ExploreData::ScoredPaths(_) => "ScoredPaths",
        }
    }
}

// ─── Filter types ───────────────────────────────────────────────────────────

/// Comparison operator for filter expressions.
#[derive(Debug, PartialEq, Eq)]
enum FilterOp {
    Eq,
    Ne,
    Gt,
    Lt,
    Ge,
    Le,
    Contains,
}

/// A single filter condition: `field op value`.
#[derive(Debug)]
struct FilterExpr {
    field: String,
    op: FilterOp,
    value: String,
}

/// Parse a filter expression from tokenised command parts.
///
/// Expected format: `["filter", "<field>", "<op>", "<value>", ...]`.
/// The value may span multiple tokens (joined with spaces).
///
/// Returns `None` if the input is too short or the operator is unrecognised.
fn parse_filter(parts: &[&str]) -> Option<FilterExpr> {
    if parts.len() < 4 {
        return None;
    }
    let field = parts[1].to_string();
    let op = match parts[2] {
        "=" | "==" => FilterOp::Eq,
        "!=" => FilterOp::Ne,
        ">" => FilterOp::Gt,
        "<" => FilterOp::Lt,
        ">=" => FilterOp::Ge,
        "<=" => FilterOp::Le,
        "contains" => FilterOp::Contains,
        _ => return None,
    };
    let value = parts[3..].join(" ");
    Some(FilterExpr { field, op, value })
}

// ─── Field extraction ───────────────────────────────────────────────────────

/// Extract a named field from a Molecule as a string.
fn molecule_field(mol: &Molecule, field: &str) -> Option<String> {
    match field {
        "id" => Some(mol.id.to_string()),
        "umi_fwd" => Some(String::from_utf8_lossy(&mol.umi_fwd).into_owned()),
        "umi_rev" => Some(String::from_utf8_lossy(&mol.umi_rev).into_owned()),
        "n_reads_fwd" => {
            let n: u8 = mol.consensus_fwd.as_ref().map_or(0, |c| c.family_size.0);
            Some(n.to_string())
        }
        "n_reads_rev" => {
            let n: u8 = mol.consensus_rev.as_ref().map_or(0, |c| c.family_size.1);
            Some(n.to_string())
        }
        "is_duplex" => {
            let fwd = mol.consensus_fwd.is_some();
            let rev = mol.consensus_rev.is_some();
            Some((fwd && rev).to_string())
        }
        "family_type" => {
            let fwd_reads = mol.consensus_fwd.as_ref().map_or(0, |c| c.family_size.0);
            let rev_reads = mol.consensus_rev.as_ref().map_or(0, |c| c.family_size.1);
            let ft = FamilyType::from_family_size((fwd_reads, rev_reads));
            Some(format!("{ft:?}"))
        }
        "seq_len" => {
            // Length of the longest available consensus.
            let len = mol
                .duplex_consensus
                .as_ref()
                .map(|c| c.sequence.len())
                .or_else(|| mol.consensus_fwd.as_ref().map(|c| c.sequence.len()))
                .or_else(|| mol.consensus_rev.as_ref().map(|c| c.sequence.len()))
                .unwrap_or(0);
            Some(len.to_string())
        }
        _ => None,
    }
}

/// Extract a named field from a KmerEntry as a string.
fn kmer_field(entry: &KmerEntry, field: &str) -> Option<String> {
    match field {
        "kmer" => Some(entry.kmer.to_string()),
        "n_molecules" => Some(entry.n_molecules.to_string()),
        "n_duplex" => Some(entry.n_duplex.to_string()),
        "n_simplex_fwd" => Some(entry.n_simplex_fwd.to_string()),
        "n_simplex_rev" => Some(entry.n_simplex_rev.to_string()),
        "min_base_error_prob" => Some(format!("{:.6}", entry.min_base_error_prob)),
        "mean_base_error_prob" => Some(format!("{:.6}", entry.mean_base_error_prob)),
        _ => None,
    }
}

/// Extract a named field from a ScoredPathRecord as a string.
fn path_field(rec: &ScoredPathRecord, field: &str) -> Option<String> {
    match field {
        "target_id" => Some(rec.target_id.clone()),
        "is_reference" => Some(rec.is_reference.to_string()),
        "min_molecules" => Some(rec.min_molecules.to_string()),
        "mean_molecules" => Some(format!("{:.2}", rec.mean_molecules)),
        "min_duplex" => Some(rec.min_duplex.to_string()),
        "mean_duplex" => Some(format!("{:.2}", rec.mean_duplex)),
        "min_variant_specific_duplex" => Some(rec.min_variant_specific_duplex.to_string()),
        "mean_variant_specific_molecules" => {
            Some(format!("{:.2}", rec.mean_variant_specific_molecules))
        }
        "min_simplex_fwd" => Some(rec.min_simplex_fwd.to_string()),
        "min_simplex_rev" => Some(rec.min_simplex_rev.to_string()),
        "mean_error_prob" => Some(format!("{:.6}", rec.mean_error_prob)),
        "seq_len" => Some(rec.sequence.len().to_string()),
        _ => None,
    }
}

/// Extract a field value from a record at a given index.
fn extract_field(data: &ExploreData, idx: usize, field: &str) -> Option<String> {
    match data {
        ExploreData::Molecules(v) => v.get(idx).and_then(|m| molecule_field(m, field)),
        ExploreData::KmerIndex(v) => v.get(idx).and_then(|e| kmer_field(e, field)),
        ExploreData::ScoredPaths(v) => v.get(idx).and_then(|r| path_field(r, field)),
    }
}

// ─── Filter application ────────────────────────────────────────────────────

/// Compare two strings using the given operator.
///
/// For numeric operators (>, <, >=, <=), both values are parsed as f64.
/// If parsing fails, the comparison returns false.
fn compare(lhs: &str, op: &FilterOp, rhs: &str) -> bool {
    match op {
        FilterOp::Eq => lhs == rhs,
        FilterOp::Ne => lhs != rhs,
        FilterOp::Contains => lhs.contains(rhs),
        FilterOp::Gt | FilterOp::Lt | FilterOp::Ge | FilterOp::Le => {
            let a: f64 = match lhs.parse() {
                Ok(v) => v,
                Err(_) => return false,
            };
            let b: f64 = match rhs.parse() {
                Ok(v) => v,
                Err(_) => return false,
            };
            match op {
                FilterOp::Gt => a > b,
                FilterOp::Lt => a < b,
                FilterOp::Ge => a >= b,
                FilterOp::Le => a <= b,
                _ => unreachable!(),
            }
        }
    }
}

/// Return indices of records matching the filter expression.
fn apply_filter(data: &ExploreData, expr: &FilterExpr) -> Vec<usize> {
    let n = data.len();
    let mut matches = Vec::new();
    for i in 0..n {
        if let Some(val) = extract_field(data, i, &expr.field) {
            if compare(&val, &expr.op, &expr.value) {
                matches.push(i);
            }
        }
    }
    matches
}

// ─── File loading ───────────────────────────────────────────────────────────

/// Load a kam bincode file and return the typed data.
fn load_file(path: &Path) -> Result<ExploreData, Box<dyn std::error::Error>> {
    let header = read_header(path)?;
    match header.file_type {
        FileType::Molecules => {
            let (_, records): (_, Vec<Molecule>) = read_bincode(path)?;
            Ok(ExploreData::Molecules(records))
        }
        FileType::KmerIndex => {
            let (_, records): (_, Vec<KmerEntry>) = read_bincode(path)?;
            Ok(ExploreData::KmerIndex(records))
        }
        FileType::ScoredPaths => {
            let (_, records): (_, Vec<ScoredPathRecord>) = read_bincode(path)?;
            Ok(ExploreData::ScoredPaths(records))
        }
        FileType::VariantCalls => Err(
            "VariantCalls files are not yet supported by explore (no bincode serialisation for VariantCall)"
                .into(),
        ),
    }
}

// ─── Table builders ─────────────────────────────────────────────────────────

/// Column headers for each data type.
fn column_headers(data: &ExploreData) -> Vec<&'static str> {
    match data {
        ExploreData::Molecules(_) => vec![
            "#",
            "id",
            "umi_fwd",
            "umi_rev",
            "family_type",
            "n_reads_fwd",
            "n_reads_rev",
            "seq_len",
        ],
        ExploreData::KmerIndex(_) => vec![
            "#",
            "kmer",
            "n_molecules",
            "n_duplex",
            "n_simplex_fwd",
            "n_simplex_rev",
            "min_base_error_prob",
            "mean_base_error_prob",
        ],
        ExploreData::ScoredPaths(_) => vec![
            "#",
            "target_id",
            "is_reference",
            "min_molecules",
            "mean_molecules",
            "min_duplex",
            "mean_duplex",
            "seq_len",
        ],
    }
}

/// Fields to extract for table rows (must match column_headers minus the "#" column).
fn row_fields(data: &ExploreData) -> Vec<&'static str> {
    match data {
        ExploreData::Molecules(_) => vec![
            "id",
            "umi_fwd",
            "umi_rev",
            "family_type",
            "n_reads_fwd",
            "n_reads_rev",
            "seq_len",
        ],
        ExploreData::KmerIndex(_) => vec![
            "kmer",
            "n_molecules",
            "n_duplex",
            "n_simplex_fwd",
            "n_simplex_rev",
            "min_base_error_prob",
            "mean_base_error_prob",
        ],
        ExploreData::ScoredPaths(_) => vec![
            "target_id",
            "is_reference",
            "min_molecules",
            "mean_molecules",
            "min_duplex",
            "mean_duplex",
            "seq_len",
        ],
    }
}

/// Build a comfy-table from the given record indices.
fn build_table(data: &ExploreData, indices: &[usize]) -> Table {
    let mut table = Table::new();
    let headers: Vec<Cell> = column_headers(data).into_iter().map(Cell::new).collect();
    table.set_header(headers);

    let fields = row_fields(data);
    for &idx in indices {
        let mut row: Vec<Cell> = vec![Cell::new(idx)];
        for f in &fields {
            let val = extract_field(data, idx, f).unwrap_or_default();
            row.push(Cell::new(val));
        }
        table.add_row(row);
    }
    table
}

// ─── Commands ───────────────────────────────────────────────────────────────

/// Print a one-line summary of the loaded data.
fn print_summary(data: &ExploreData, file_size: u64) {
    println!(
        "File type: {}  |  Records: {}  |  File size: {}",
        data.type_label(),
        data.len(),
        format_size(file_size),
    );
}

/// Format a byte count as a human-readable string.
fn format_size(bytes: u64) -> String {
    const KIB: u64 = 1024;
    const MIB: u64 = KIB * 1024;
    const GIB: u64 = MIB * 1024;
    if bytes >= GIB {
        format!("{:.2} GiB", bytes as f64 / GIB as f64)
    } else if bytes >= MIB {
        format!("{:.2} MiB", bytes as f64 / MIB as f64)
    } else if bytes >= KIB {
        format!("{:.2} KiB", bytes as f64 / KIB as f64)
    } else {
        format!("{bytes} B")
    }
}

/// Show the first N records in a table.
fn cmd_head(data: &ExploreData, parts: &[&str]) {
    let n: usize = parts
        .get(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(10)
        .min(data.len());
    let indices: Vec<usize> = (0..n).collect();
    println!("{}", build_table(data, &indices));
}

/// Show a single record in full detail.
fn cmd_show(data: &ExploreData, parts: &[&str]) {
    let idx: usize = match parts.get(1).and_then(|s| s.parse().ok()) {
        Some(i) => i,
        None => {
            eprintln!("Usage: show <index>");
            return;
        }
    };
    if idx >= data.len() {
        eprintln!("Index {idx} out of range (0..{})", data.len());
        return;
    }

    match data {
        ExploreData::Molecules(v) => show_molecule(&v[idx], idx),
        ExploreData::KmerIndex(v) => show_kmer_entry(&v[idx], idx),
        ExploreData::ScoredPaths(v) => show_scored_path(&v[idx], idx),
    }
}

/// Print all fields of a single Molecule.
fn show_molecule(mol: &Molecule, idx: usize) {
    let fwd_reads = mol.consensus_fwd.as_ref().map_or(0, |c| c.family_size.0);
    let rev_reads = mol.consensus_rev.as_ref().map_or(0, |c| c.family_size.1);
    let ft = FamilyType::from_family_size((fwd_reads, rev_reads));

    println!("── Molecule #{idx} ──");
    println!("  id:           {}", mol.id);
    println!(
        "  umi_fwd:      {}",
        String::from_utf8_lossy(&mol.umi_fwd)
    );
    println!(
        "  umi_rev:      {}",
        String::from_utf8_lossy(&mol.umi_rev)
    );
    println!("  family_type:  {ft:?}");
    println!("  n_reads_fwd:  {fwd_reads}");
    println!("  n_reads_rev:  {rev_reads}");

    if let Some(ref c) = mol.consensus_fwd {
        println!("  consensus_fwd: {} bp", c.sequence.len());
    } else {
        println!("  consensus_fwd: (none)");
    }
    if let Some(ref c) = mol.consensus_rev {
        println!("  consensus_rev: {} bp", c.sequence.len());
    } else {
        println!("  consensus_rev: (none)");
    }
    if let Some(ref c) = mol.duplex_consensus {
        println!("  duplex_consensus: {} bp", c.sequence.len());
    } else {
        println!("  duplex_consensus: (none)");
    }
    println!(
        "  evidence:     {}",
        if mol.evidence.is_some() {
            "present"
        } else {
            "(none)"
        }
    );
}

/// Print all fields of a single KmerEntry.
fn show_kmer_entry(entry: &KmerEntry, idx: usize) {
    println!("── KmerEntry #{idx} ──");
    println!("  kmer (encoded): {}", entry.kmer);
    println!("  n_molecules:    {}", entry.n_molecules);
    println!("  n_duplex:       {}", entry.n_duplex);
    println!("  n_simplex_fwd:  {}", entry.n_simplex_fwd);
    println!("  n_simplex_rev:  {}", entry.n_simplex_rev);
    println!("  min_base_error_prob:  {:.6}", entry.min_base_error_prob);
    println!("  mean_base_error_prob: {:.6}", entry.mean_base_error_prob);
}

/// Print all fields of a single ScoredPathRecord.
fn show_scored_path(rec: &ScoredPathRecord, idx: usize) {
    println!("── ScoredPathRecord #{idx} ──");
    println!("  target_id:     {}", rec.target_id);
    println!("  is_reference:  {}", rec.is_reference);
    println!("  sequence:      {} bp", rec.sequence.len());
    println!("  min_molecules: {}", rec.min_molecules);
    println!("  mean_molecules: {:.2}", rec.mean_molecules);
    println!("  min_duplex:    {}", rec.min_duplex);
    println!("  mean_duplex:   {:.2}", rec.mean_duplex);
    println!(
        "  min_variant_specific_duplex:      {}",
        rec.min_variant_specific_duplex
    );
    println!(
        "  mean_variant_specific_molecules:  {:.2}",
        rec.mean_variant_specific_molecules
    );
    println!("  min_simplex_fwd: {}", rec.min_simplex_fwd);
    println!("  min_simplex_rev: {}", rec.min_simplex_rev);
    println!("  mean_error_prob: {:.6}", rec.mean_error_prob);
}

/// Run the filter command: parse the expression, apply it, and show results.
fn cmd_filter(data: &ExploreData, parts: &[&str]) {
    let expr = match parse_filter(parts) {
        Some(e) => e,
        None => {
            eprintln!("Usage: filter <field> <op> <value>");
            eprintln!("  Operators: =, !=, >, <, >=, <=, contains");
            eprintln!("  Example:   filter n_molecules > 5");
            return;
        }
    };

    // Check that the field name is valid for this data type.
    if extract_field(data, 0, &expr.field).is_none() && data.len() > 0 {
        eprintln!(
            "Unknown field '{}' for {} records.",
            expr.field,
            data.type_label()
        );
        eprintln!("Available fields: {:?}", row_fields(data));
        return;
    }

    let matches = apply_filter(data, &expr);
    println!(
        "{} of {} records match: {} {:?} {}",
        matches.len(),
        data.len(),
        expr.field,
        expr.op,
        expr.value,
    );

    if !matches.is_empty() {
        let show: Vec<usize> = matches.iter().copied().take(10).collect();
        println!("{}", build_table(data, &show));
        if matches.len() > 10 {
            println!("  ... and {} more", matches.len() - 10);
        }
    }
}

/// Print summary statistics for the loaded data.
fn cmd_stats(data: &ExploreData) {
    println!("── Statistics ──");
    println!("  Records: {}", data.len());

    match data {
        ExploreData::Molecules(v) => stats_molecules(v),
        ExploreData::KmerIndex(v) => stats_kmer_index(v),
        ExploreData::ScoredPaths(v) => stats_scored_paths(v),
    }
}

/// Molecule-specific statistics.
fn stats_molecules(mols: &[Molecule]) {
    if mols.is_empty() {
        return;
    }

    let mut n_duplex: usize = 0;
    let mut n_simplex_fwd: usize = 0;
    let mut n_simplex_rev: usize = 0;
    let mut n_singleton: usize = 0;
    let mut total_reads: u64 = 0;

    for mol in mols {
        let fwd = mol.consensus_fwd.as_ref().map_or(0, |c| c.family_size.0);
        let rev = mol.consensus_rev.as_ref().map_or(0, |c| c.family_size.1);
        total_reads += u64::from(fwd) + u64::from(rev);
        match FamilyType::from_family_size((fwd, rev)) {
            FamilyType::Duplex => n_duplex += 1,
            FamilyType::SimplexFwd => n_simplex_fwd += 1,
            FamilyType::SimplexRev => n_simplex_rev += 1,
            FamilyType::Singleton => n_singleton += 1,
        }
    }

    println!("  Duplex:      {n_duplex}");
    println!("  SimplexFwd:  {n_simplex_fwd}");
    println!("  SimplexRev:  {n_simplex_rev}");
    println!("  Singleton:   {n_singleton}");
    println!("  Total reads: {total_reads}");
    println!(
        "  Mean reads/molecule: {:.1}",
        total_reads as f64 / mols.len() as f64
    );
}

/// KmerIndex-specific statistics.
fn stats_kmer_index(entries: &[KmerEntry]) {
    if entries.is_empty() {
        return;
    }

    let mut total_mols: u64 = 0;
    let mut total_duplex: u64 = 0;
    let mut max_mols: u32 = 0;
    let mut min_mols: u32 = u32::MAX;

    for e in entries {
        total_mols += u64::from(e.n_molecules);
        total_duplex += u64::from(e.n_duplex);
        max_mols = max_mols.max(e.n_molecules);
        min_mols = min_mols.min(e.n_molecules);
    }

    println!("  Distinct k-mers: {}", entries.len());
    println!(
        "  Molecules per k-mer: min={min_mols}, max={max_mols}, mean={:.1}",
        total_mols as f64 / entries.len() as f64
    );
    println!(
        "  Duplex fraction: {:.1}%",
        if total_mols > 0 {
            total_duplex as f64 / total_mols as f64 * 100.0
        } else {
            0.0
        }
    );
}

/// ScoredPaths-specific statistics.
fn stats_scored_paths(recs: &[ScoredPathRecord]) {
    if recs.is_empty() {
        return;
    }

    let n_ref = recs.iter().filter(|r| r.is_reference).count();
    let n_alt = recs.len() - n_ref;

    let unique_targets: std::collections::HashSet<&str> =
        recs.iter().map(|r| r.target_id.as_str()).collect();

    println!("  Unique targets: {}", unique_targets.len());
    println!("  Reference paths: {n_ref}");
    println!("  Alternate paths: {n_alt}");

    // Molecule support distribution across alt paths.
    let alt_mols: Vec<u32> = recs
        .iter()
        .filter(|r| !r.is_reference)
        .map(|r| r.min_molecules)
        .collect();
    if !alt_mols.is_empty() {
        let min = alt_mols.iter().copied().min().unwrap_or(0);
        let max = alt_mols.iter().copied().max().unwrap_or(0);
        let mean = alt_mols.iter().map(|&v| f64::from(v)).sum::<f64>() / alt_mols.len() as f64;
        println!("  Alt min_molecules: min={min}, max={max}, mean={mean:.1}");
    }
}

/// Export all records (or a filtered subset) to a TSV file.
fn cmd_export(data: &ExploreData, parts: &[&str]) {
    let path = match parts.get(1) {
        Some(p) => p,
        None => {
            eprintln!("Usage: export <path.tsv>");
            return;
        }
    };

    let fields = row_fields(data);
    let n = data.len();

    let result = (|| -> Result<usize, Box<dyn std::error::Error>> {
        let file = fs::File::create(path)?;
        let mut w = std::io::BufWriter::new(file);

        // Header row.
        writeln!(w, "{}", fields.join("\t"))?;

        // Data rows.
        for i in 0..n {
            let vals: Vec<String> = fields
                .iter()
                .map(|f| extract_field(data, i, f).unwrap_or_default())
                .collect();
            writeln!(w, "{}", vals.join("\t"))?;
        }

        w.flush()?;
        Ok(n)
    })();

    match result {
        Ok(count) => println!("Exported {count} records to {path}"),
        Err(e) => eprintln!("Export failed: {e}"),
    }
}

/// Print the help text listing all available commands.
fn print_help() {
    println!("Commands:");
    println!("  summary             Show file type, record count, and file size");
    println!("  head [N]            Show the first N records (default 10)");
    println!("  show <index>        Show one record in full detail");
    println!("  filter <field> <op> <value>");
    println!("                      Filter records by a condition and show matches");
    println!("                      Operators: =, !=, >, <, >=, <=, contains");
    println!("  stats               Summary statistics for the loaded data");
    println!("  export <path.tsv>   Export all records to a TSV file");
    println!("  help                Show this help");
    println!("  quit / exit / q     Exit the REPL");
}

// ─── Entry point ────────────────────────────────────────────────────────────

/// Run the `kam explore` subcommand.
///
/// Loads the specified bincode file, prints a summary, then enters an
/// interactive REPL where users can query and inspect the data.
///
/// # Errors
///
/// Returns an error if the file cannot be opened or deserialised.
pub fn run_explore(args: ExploreArgs) -> Result<(), Box<dyn std::error::Error>> {
    let file_size = fs::metadata(&args.file)
        .map(|m| m.len())
        .unwrap_or(0);

    let data = load_file(&args.file)?;
    print_summary(&data, file_size);
    println!("Type 'help' for available commands.\n");

    let mut rl = rustyline::DefaultEditor::new()?;
    loop {
        let line = match rl.readline("kam> ") {
            Ok(line) => line,
            Err(ReadlineError::Interrupted | ReadlineError::Eof) => break,
            Err(e) => return Err(e.into()),
        };
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let _ = rl.add_history_entry(trimmed);

        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        match parts[0] {
            "summary" => print_summary(&data, file_size),
            "head" => cmd_head(&data, &parts),
            "show" => cmd_show(&data, &parts),
            "filter" => cmd_filter(&data, &parts),
            "stats" => cmd_stats(&data),
            "export" => cmd_export(&data, &parts),
            "help" => print_help(),
            "quit" | "exit" | "q" => break,
            cmd => eprintln!("Unknown command: {cmd}. Type 'help' for commands."),
        }
    }

    Ok(())
}

// ─── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── parse_filter ────────────────────────────────────────────────────────

    #[test]
    fn parse_filter_eq() {
        let parts = vec!["filter", "variant_type", "=", "SNV"];
        let expr = parse_filter(&parts).expect("should parse");
        assert_eq!(expr.field, "variant_type");
        assert!(matches!(expr.op, FilterOp::Eq));
        assert_eq!(expr.value, "SNV");
    }

    #[test]
    fn parse_filter_double_eq() {
        let parts = vec!["filter", "variant_type", "==", "SNV"];
        let expr = parse_filter(&parts).expect("should parse");
        assert!(matches!(expr.op, FilterOp::Eq));
    }

    #[test]
    fn parse_filter_ne() {
        let parts = vec!["filter", "family_type", "!=", "Singleton"];
        let expr = parse_filter(&parts).expect("should parse");
        assert!(matches!(expr.op, FilterOp::Ne));
        assert_eq!(expr.value, "Singleton");
    }

    #[test]
    fn parse_filter_gt() {
        let parts = vec!["filter", "vaf", ">", "0.01"];
        let expr = parse_filter(&parts).expect("should parse");
        assert_eq!(expr.field, "vaf");
        assert!(matches!(expr.op, FilterOp::Gt));
        assert_eq!(expr.value, "0.01");
    }

    #[test]
    fn parse_filter_lt() {
        let parts = vec!["filter", "n_molecules", "<", "10"];
        let expr = parse_filter(&parts).expect("should parse");
        assert!(matches!(expr.op, FilterOp::Lt));
    }

    #[test]
    fn parse_filter_ge() {
        let parts = vec!["filter", "n_duplex", ">=", "2"];
        let expr = parse_filter(&parts).expect("should parse");
        assert!(matches!(expr.op, FilterOp::Ge));
    }

    #[test]
    fn parse_filter_le() {
        let parts = vec!["filter", "mean_error_prob", "<=", "0.001"];
        let expr = parse_filter(&parts).expect("should parse");
        assert!(matches!(expr.op, FilterOp::Le));
    }

    #[test]
    fn parse_filter_contains() {
        let parts = vec!["filter", "target_id", "contains", "BRCA"];
        let expr = parse_filter(&parts).expect("should parse");
        assert!(matches!(expr.op, FilterOp::Contains));
        assert_eq!(expr.value, "BRCA");
    }

    #[test]
    fn parse_filter_value_with_spaces() {
        let parts = vec!["filter", "target_id", "contains", "BRCA", "exon", "3"];
        let expr = parse_filter(&parts).expect("should parse");
        assert_eq!(expr.value, "BRCA exon 3");
    }

    #[test]
    fn parse_filter_too_short_returns_none() {
        let parts = vec!["filter", "vaf"];
        assert!(parse_filter(&parts).is_none());
    }

    #[test]
    fn parse_filter_unknown_op_returns_none() {
        let parts = vec!["filter", "vaf", "~", "0.01"];
        assert!(parse_filter(&parts).is_none());
    }

    // ── compare ─────────────────────────────────────────────────────────────

    #[test]
    fn compare_eq_string() {
        assert!(compare("SNV", &FilterOp::Eq, "SNV"));
        assert!(!compare("SNV", &FilterOp::Eq, "MNV"));
    }

    #[test]
    fn compare_ne_string() {
        assert!(compare("SNV", &FilterOp::Ne, "MNV"));
        assert!(!compare("SNV", &FilterOp::Ne, "SNV"));
    }

    #[test]
    fn compare_gt_numeric() {
        assert!(compare("0.05", &FilterOp::Gt, "0.01"));
        assert!(!compare("0.01", &FilterOp::Gt, "0.05"));
    }

    #[test]
    fn compare_lt_numeric() {
        assert!(compare("3", &FilterOp::Lt, "10"));
        assert!(!compare("10", &FilterOp::Lt, "3"));
    }

    #[test]
    fn compare_ge_numeric() {
        assert!(compare("5", &FilterOp::Ge, "5"));
        assert!(compare("6", &FilterOp::Ge, "5"));
        assert!(!compare("4", &FilterOp::Ge, "5"));
    }

    #[test]
    fn compare_le_numeric() {
        assert!(compare("5", &FilterOp::Le, "5"));
        assert!(compare("4", &FilterOp::Le, "5"));
        assert!(!compare("6", &FilterOp::Le, "5"));
    }

    #[test]
    fn compare_contains() {
        assert!(compare("BRCA1_exon10", &FilterOp::Contains, "BRCA"));
        assert!(!compare("TP53_exon5", &FilterOp::Contains, "BRCA"));
    }

    #[test]
    fn compare_numeric_op_on_non_numeric_returns_false() {
        assert!(!compare("abc", &FilterOp::Gt, "5"));
        assert!(!compare("5", &FilterOp::Lt, "xyz"));
    }

    // ── format_size ─────────────────────────────────────────────────────────

    #[test]
    fn format_size_bytes() {
        assert_eq!(format_size(512), "512 B");
    }

    #[test]
    fn format_size_kib() {
        assert_eq!(format_size(2048), "2.00 KiB");
    }

    #[test]
    fn format_size_mib() {
        assert_eq!(format_size(5 * 1024 * 1024), "5.00 MiB");
    }

    #[test]
    fn format_size_gib() {
        assert_eq!(format_size(3 * 1024 * 1024 * 1024), "3.00 GiB");
    }

    // ── apply_filter on synthetic data ──────────────────────────────────────

    fn make_test_kmer_entries() -> ExploreData {
        ExploreData::KmerIndex(vec![
            KmerEntry {
                kmer: 1,
                n_molecules: 10,
                n_duplex: 5,
                n_simplex_fwd: 3,
                n_simplex_rev: 2,
                min_base_error_prob: 0.001,
                mean_base_error_prob: 0.002,
            },
            KmerEntry {
                kmer: 2,
                n_molecules: 3,
                n_duplex: 1,
                n_simplex_fwd: 1,
                n_simplex_rev: 1,
                min_base_error_prob: 0.005,
                mean_base_error_prob: 0.01,
            },
            KmerEntry {
                kmer: 3,
                n_molecules: 20,
                n_duplex: 12,
                n_simplex_fwd: 4,
                n_simplex_rev: 4,
                min_base_error_prob: 0.0001,
                mean_base_error_prob: 0.0005,
            },
        ])
    }

    #[test]
    fn apply_filter_gt_on_kmer_entries() {
        let data = make_test_kmer_entries();
        let expr = FilterExpr {
            field: "n_molecules".to_string(),
            op: FilterOp::Gt,
            value: "5".to_string(),
        };
        let matches = apply_filter(&data, &expr);
        assert_eq!(matches, vec![0, 2]);
    }

    #[test]
    fn apply_filter_eq_on_kmer_entries() {
        let data = make_test_kmer_entries();
        let expr = FilterExpr {
            field: "n_duplex".to_string(),
            op: FilterOp::Eq,
            value: "1".to_string(),
        };
        let matches = apply_filter(&data, &expr);
        assert_eq!(matches, vec![1]);
    }

    #[test]
    fn apply_filter_unknown_field_returns_empty() {
        let data = make_test_kmer_entries();
        let expr = FilterExpr {
            field: "nonexistent".to_string(),
            op: FilterOp::Eq,
            value: "foo".to_string(),
        };
        let matches = apply_filter(&data, &expr);
        assert!(matches.is_empty());
    }

    // ── Molecule field extraction ───────────────────────────────────────────

    fn make_test_molecule() -> Molecule {
        use kam_core::molecule::ConsensusRead;
        Molecule {
            id: 42,
            umi_fwd: b"ACGTA".to_vec(),
            umi_rev: b"TGCAT".to_vec(),
            consensus_fwd: Some(ConsensusRead {
                sequence: b"ACGTACGT".to_vec(),
                per_base_error_prob: vec![0.001; 8],
                per_base_strand_support: vec![(2, 0); 8],
                family_size: (3, 0),
            }),
            consensus_rev: Some(ConsensusRead {
                sequence: b"ACGTACGT".to_vec(),
                per_base_error_prob: vec![0.001; 8],
                per_base_strand_support: vec![(0, 2); 8],
                family_size: (0, 2),
            }),
            duplex_consensus: None,
            evidence: None,
        }
    }

    #[test]
    fn molecule_field_id() {
        let mol = make_test_molecule();
        assert_eq!(molecule_field(&mol, "id"), Some("42".to_string()));
    }

    #[test]
    fn molecule_field_umi_fwd() {
        let mol = make_test_molecule();
        assert_eq!(molecule_field(&mol, "umi_fwd"), Some("ACGTA".to_string()));
    }

    #[test]
    fn molecule_field_family_type_duplex() {
        let mol = make_test_molecule();
        assert_eq!(
            molecule_field(&mol, "family_type"),
            Some("Duplex".to_string())
        );
    }

    #[test]
    fn molecule_field_unknown_returns_none() {
        let mol = make_test_molecule();
        assert_eq!(molecule_field(&mol, "nonexistent"), None);
    }
}
