use std::fmt::Write;
use std::fs;
use std::path::Path;

/// Resource and timing metrics for a single pipeline stage.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StageMetrics {
    pub stage: &'static str,
    pub elapsed_ms: u64,
    pub peak_rss_mb: u64,
    pub cpu_time_ms: u64,
}

/// Aggregate metrics across all pipeline stages in a run.
#[derive(Debug, Clone, serde::Serialize)]
pub struct RunMetrics {
    pub stages: Vec<StageMetrics>,
    pub peak_rss_mb: u64,
    pub total_cpu_time_ms: u64,
    pub total_elapsed_ms: u64,
}

/// Read peak virtual memory from /proc/self/status.
///
/// Parses the `VmPeak:` line and returns the value in megabytes.
/// Returns 0 if /proc is not available or the line cannot be parsed.
pub fn read_peak_rss_mb() -> u64 {
    let status = match fs::read_to_string("/proc/self/status") {
        Ok(s) => s,
        Err(_) => return 0,
    };

    for line in status.lines() {
        if let Some(val) = line.strip_prefix("VmPeak:") {
            // Expected format: "VmPeak:   12345 kB"
            let val = val.trim();
            if let Some(kb_str) = val.strip_suffix(" kB") {
                if let Ok(kb) = kb_str.trim().parse::<u64>() {
                    return kb / 1024;
                }
            }
        }
    }

    0
}

/// Read CPU time (user + system) from /proc/self/stat.
///
/// Fields 14 (utime) and 15 (stime) are in clock ticks.
/// Divides by 100 as an approximation of CLK_TCK for an approximate
/// millisecond value. Returns 0 if /proc is not available or the fields
/// cannot be parsed.
pub fn read_cpu_time_ms() -> u64 {
    let stat = match fs::read_to_string("/proc/self/stat") {
        Ok(s) => s,
        Err(_) => return 0,
    };

    // /proc/self/stat has a tricky format: the second field (comm) is in
    // parentheses and can contain spaces and closing parens. We need to find
    // the last ')' and skip the space after it to get to field 3.
    let after_comm = match stat.rfind(") ") {
        Some(pos) => &stat[pos + 2..],
        None => return 0,
    };

    let fields: Vec<&str> = after_comm.split_whitespace().collect();

    // Fields 14 and 15 are zero-indexed at indices 12 and 13.
    let utime = fields
        .get(12)
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(0);
    let stime = fields
        .get(13)
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(0);

    // Convert clock ticks to milliseconds (assuming 100 Hz).
    // Use saturating arithmetic in case of overflow.
    (utime.saturating_add(stime)).saturating_mul(10)
}

/// A simple scope-guarded timer that records stage metrics.
///
/// Create one at the start of a pipeline stage and call [`finish`](StageTimer::finish)
/// at the end to produce [`StageMetrics`].
pub struct StageTimer {
    stage: &'static str,
    start: std::time::Instant,
}

impl StageTimer {
    /// Start timing a new stage.
    pub fn new(stage: &'static str) -> Self {
        Self {
            stage,
            start: std::time::Instant::now(),
        }
    }

    /// Stop the timer and return metrics for this stage.
    ///
    /// The snapshot of RSS and CPU time is taken at the moment of finishing.
    pub fn finish(&mut self) -> StageMetrics {
        let elapsed_ms = self.start.elapsed().as_millis() as u64;
        StageMetrics {
            stage: self.stage,
            elapsed_ms,
            peak_rss_mb: read_peak_rss_mb(),
            cpu_time_ms: read_cpu_time_ms(),
        }
    }
}

/// Format stage metrics as a human-readable table.
///
/// The table includes each stage's elapsed time, CPU time, peak RSS, and a
/// summary row with totals.
pub fn format_metrics_table(
    stages: &[StageMetrics],
    total_elapsed_ms: u64,
    total_cpu_ms: u64,
    peak_rss_mb: u64,
) -> String {
    // Stage column width: at least 5 ("Stage"), or the longest stage name.
    let stage_w = stages
        .iter()
        .map(|s| s.stage.len())
        .max()
        .unwrap_or(0)
        .max(5);

    let num_w = 10; // width for numeric columns

    let mut table = String::new();

    // Header row.
    let hdr_stage = pad_right("Stage", stage_w);
    writeln!(
        table,
        "{hdr_stage}  {:>num_w$}  {:>num_w$}  {:>num_w$}",
        "Elapsed", "CPU", "Peak RSS",
    )
    .ok();

    // Divider.
    let divider_len = stage_w + 2 + num_w + 2 + num_w + 2 + num_w;
    writeln!(table, "{}", "-".repeat(divider_len)).ok();

    // Stage rows.
    for s in stages {
        let name = pad_right(s.stage, stage_w);
        writeln!(
            table,
            "{name}  {:>num_w$}  {:>num_w$}  {:>num_w$}",
            format_ms(s.elapsed_ms),
            format_ms(s.cpu_time_ms),
            format_mb(s.peak_rss_mb),
        )
        .ok();
    }

    // Totals separator.
    writeln!(table, "{}", "-".repeat(divider_len)).ok();

    // Totals row.
    let total_name = pad_right("Total", stage_w);
    writeln!(
        table,
        "{total_name}  {:>num_w$}  {:>num_w$}  {:>num_w$}",
        format_ms(total_elapsed_ms),
        format_ms(total_cpu_ms),
        format_mb(peak_rss_mb),
    )
    .ok();

    table
}

/// Format a millisecond value with "ms" suffix.
fn format_ms(ms: u64) -> String {
    format!("{}ms", ms)
}

/// Format an RSS value with "MB" suffix.
fn format_mb(mb: u64) -> String {
    format!("{}MB", mb)
}

/// Write metrics to a JSON file at the given path.
///
/// # Errors
///
/// Returns an error if the file cannot be written or serialisation fails.
pub fn write_metrics_json(
    metrics: &RunMetrics,
    path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let json = serde_json::to_string_pretty(metrics)?;
    fs::write(path, json)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Pad a string to the given width with spaces on the right.
fn pad_right(s: &str, width: usize) -> String {
    if s.len() >= width {
        s.to_string()
    } else {
        let mut result = String::with_capacity(width);
        result.push_str(s);
        for _ in 0..(width - s.len()) {
            result.push(' ');
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // read_peak_rss_mb
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_peak_rss_does_not_crash() {
        let _ = read_peak_rss_mb();
    }

    // -----------------------------------------------------------------------
    // read_cpu_time_ms
    // -----------------------------------------------------------------------

    #[test]
    fn test_read_cpu_time_does_not_crash() {
        let _ = read_cpu_time_ms();
    }

    // -----------------------------------------------------------------------
    // StageTimer
    // -----------------------------------------------------------------------

    #[test]
    fn test_stage_timer_finish_returns_metrics() {
        let mut timer = StageTimer::new("test_stage");
        // Burn a little CPU to ensure non-zero elapsed time.
        let mut x: u64 = 0;
        for i in 0..1_000_000 {
            x = x.wrapping_add(i);
        }
        std::hint::black_box(x);

        let metrics = timer.finish();
        assert_eq!(metrics.stage, "test_stage");
        assert!(
            metrics.elapsed_ms > 0,
            "expected positive elapsed time, got {}",
            metrics.elapsed_ms
        );
    }

    #[test]
    fn test_stage_timer_multiple_stages() {
        let mut t1 = StageTimer::new("stage_a");
        let mut t2 = StageTimer::new("stage_b");

        let a = t1.finish();
        let _ = t2.finish();

        assert_eq!(a.stage, "stage_a");
    }

    // -----------------------------------------------------------------------
    // format_metrics_table
    // -----------------------------------------------------------------------

    #[test]
    fn test_format_metrics_table_empty() {
        let table = format_metrics_table(&[], 0, 0, 0);
        assert!(table.contains("Stage"));
        assert!(table.contains("Total"));
        assert!(table.contains("Elapsed"));
        assert!(table.contains("CPU"));
        assert!(table.contains("Peak RSS"));
    }

    #[test]
    fn test_format_metrics_table_single_stage() {
        let stages = vec![StageMetrics {
            stage: "index",
            elapsed_ms: 1500,
            peak_rss_mb: 64,
            cpu_time_ms: 1200,
        }];
        let table = format_metrics_table(&stages, 1500, 1200, 64);
        assert!(table.contains("index"));
        assert!(table.contains("1500ms"));
        assert!(table.contains("64MB"));
        assert!(table.contains("Total"));
    }

    #[test]
    fn test_format_metrics_table_aligns_columns() {
        let stages = vec![
            StageMetrics {
                stage: "assemble",
                elapsed_ms: 54321,
                peak_rss_mb: 1024,
                cpu_time_ms: 50000,
            },
            StageMetrics {
                stage: "pathfind",
                elapsed_ms: 1234,
                peak_rss_mb: 256,
                cpu_time_ms: 1000,
            },
        ];
        let table = format_metrics_table(&stages, 55555, 51000, 1024);
        assert!(table.contains("assemble"));
        assert!(table.contains("pathfind"));
        assert!(table.contains("54321ms"));
        assert!(table.contains("1234ms"));
    }

    // -----------------------------------------------------------------------
    // write_metrics_json
    // -----------------------------------------------------------------------

    #[test]
    fn test_write_metrics_json_round_trip() {
        let metrics = RunMetrics {
            stages: vec![StageMetrics {
                stage: "test",
                elapsed_ms: 100,
                peak_rss_mb: 50,
                cpu_time_ms: 80,
            }],
            peak_rss_mb: 50,
            total_cpu_time_ms: 80,
            total_elapsed_ms: 100,
        };

        let dir = std::env::temp_dir();
        let path = dir.join("kam_test_metrics.json");
        let _ = fs::remove_file(&path);

        let result = write_metrics_json(&metrics, &path);
        assert!(result.is_ok(), "write failed: {:?}", result);

        let contents = fs::read_to_string(&path).expect("failed to read back");
        assert!(contents.contains("\"stage\": \"test\""));
        assert!(contents.contains("\"elapsed_ms\": 100"));
        assert!(contents.contains("\"total_elapsed_ms\": 100"));

        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_write_metrics_json_invalid_path() {
        let metrics = RunMetrics {
            stages: vec![],
            peak_rss_mb: 0,
            total_cpu_time_ms: 0,
            total_elapsed_ms: 0,
        };
        let result = write_metrics_json(&metrics, Path::new("/nonexistent_dir/file.json"));
        assert!(result.is_err());
    }
}
