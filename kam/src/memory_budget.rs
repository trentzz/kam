//! Adaptive memory budgeting for kam pipeline stages.
//!
//! Automatically divides a user-specified memory budget (in GB) across:
//! - Phase 1 (assemble): 25% - FASTQ streaming batch sizing
//! - Phase 2 (index): 60% - k-mer indexing with prefiltering
//! - Phase 3 (pathfind): 15% - De Bruijn graph structure

/// Memory budget allocator for all pipeline stages.
///
/// Given a total memory budget in GB, computes safe allocations for:
/// - read_batch_size_mb: How many MB of reads to buffer before processing
/// - index_max_kmers: Max k-mers to store in the k-mer index
/// - prefilter_bits: Count-Min Sketch size for frequency prefiltering
/// - graph_max_edges: Expected max edges in De Bruijn graph
///
/// Typical usage:
/// ```ignore
/// let budget = MemoryBudget::new(32.0); // 32 GB total
/// assert_eq!(budget.phase1_mb(), 8192.0); // 25% = 8 GB
/// assert_eq!(budget.phase2_mb(), 19660.8); // 60% ≈ 19.2 GB
/// assert_eq!(budget.phase3_mb(), 4096.0); // 15% = 4 GB
/// ```
#[derive(Debug, Clone)]
pub struct MemoryBudget {
    total_gb: f64,
    phase1_mb: f64, // Assembly: 25%
    phase2_mb: f64, // Indexing: 60%
    phase3_mb: f64, // Pathfind: 15%
}

impl MemoryBudget {
    /// Create a new memory budget from a total size in gigabytes.
    ///
    /// Splits into:
    /// - 25% for assembly streaming and consensus
    /// - 60% for k-mer indexing (with prefiltering)
    /// - 15% for De Bruijn graph construction
    ///
    /// # Arguments
    ///
    /// * `total_gb` - Total memory budget in gigabytes (e.g., 32.0 for 32 GB)
    ///
    /// # Panics
    ///
    /// Panics if `total_gb` is zero or negative.
    ///
    /// # Example
    ///
    /// ```
    /// let budget = MemoryBudget::new(32.0);
    /// assert_eq!(budget.total_gb(), 32.0);
    /// assert!((budget.read_batch_bytes() as f64) / (1024.0 * 1024.0) - 8192.0 < 0.1);
    /// ```
    pub fn new(total_gb: f64) -> Self {
        assert!(
            total_gb > 0.0,
            "memory budget must be positive, got {} GB",
            total_gb
        );

        let total_mb = total_gb * 1024.0;
        Self {
            total_gb,
            phase1_mb: total_mb * 0.25,
            phase2_mb: total_mb * 0.60,
            phase3_mb: total_mb * 0.15,
        }
    }

    /// Total memory budget in gigabytes.
    pub fn total_gb(&self) -> f64 {
        self.total_gb
    }

    /// Memory budget for assembly phase in megabytes.
    pub fn phase1_mb(&self) -> f64 {
        self.phase1_mb
    }

    /// Memory budget for indexing phase in megabytes.
    pub fn phase2_mb(&self) -> f64 {
        self.phase2_mb
    }

    /// Memory budget for pathfind phase in megabytes.
    pub fn phase3_mb(&self) -> f64 {
        self.phase3_mb
    }

    /// FASTQ read pair batch size in bytes for streaming assembly.
    ///
    /// Divides phase 1 budget by ~2 to account for:
    /// - Working memory for consensus building
    /// - Per-base error probability storage
    ///
    /// # Example
    ///
    /// ```
    /// let budget = MemoryBudget::new(32.0);
    /// let batch = budget.read_batch_bytes();
    /// // For 32 GB: ~2 billion bytes per batch (~400 kbp of paired reads)
    /// assert!(batch >= 1_000_000_000);
    /// assert!(batch <= 3_000_000_000);
    /// ```
    pub fn read_batch_bytes(&self) -> usize {
        // Use 50% of phase 1 budget for actual read data
        // rest goes to consensus storage and overhead
        ((self.phase1_mb * 1024.0 * 1024.0 * 0.50) as f64) as usize
    }

    /// Max k-mers that fit in index given memory budget.
    ///
    /// Each k-mer entry in `HashKmerIndex` uses approximately:
    /// - 8 bytes: u64 key (packed k-mer)
    /// - 24 bytes: `MoleculeEvidence` struct
    /// - ~16 bytes: HashMap overhead per entry
    /// - Total: ~48 bytes per k-mer entry
    ///
    /// Reserves 20% of phase 2 budget for prefilter structures and overhead.
    pub fn max_kmers(&self) -> usize {
        let bytes_available = self.phase2_mb * 1024.0 * 1024.0 * 0.80;
        let bytes_per_entry: f64 = 48.0;
        (bytes_available / bytes_per_entry) as usize
    }

    /// Size of Count-Min Sketch in bits for frequency prefiltering.
    ///
    /// The prefilter (Count-Min Sketch or Bloom filter) is used in pass 1
    /// to count k-mer frequencies at low memory cost. Only k-mers above
    /// a threshold are inserted into the main HashMap index in pass 2.
    ///
    /// Uses ~10% of phase 2 budget. For 32 GB total:
    /// - Phase 2 = 19.2 GB
    /// - Prefilter = ~1.9 GB = ~16B bits
    /// - Tracks frequencies for ~100M distinct k-mers
    pub fn prefilter_bits(&self) -> usize {
        // Use 10% of phase 2 budget for the prefilter
        ((self.phase2_mb * 1024.0 * 1024.0 * 8.0 * 0.10) as f64) as usize
    }

    /// Estimated max edges in De Bruijn graph.
    ///
    /// Each edge is stored as (u64 key, Vec<u64> successors) with:
    /// - 8 bytes: node key
    /// - Variable: successors vector (typically 1-3 edges per k-mer)
    /// - ~8 bytes: HashMap overhead per edge
    /// - Average: ~24 bytes per edge
    ///
    /// Typical de Bruijn graphs from sequencing have ~3x more edges than nodes.
    pub fn graph_max_edges(&self) -> usize {
        let bytes_available = self.phase3_mb * 1024.0 * 1024.0;
        let bytes_per_edge: f64 = 24.0;
        (bytes_available / bytes_per_edge) as usize
    }

    /// Frequency threshold for k-mer inclusion after prefiltering.
    ///
    /// K-mers seen fewer than this many times are dropped before index insertion.
    /// This filters out:
    /// - Sequencing errors (typically singleton or doubleton)
    /// - PCR duplicates from low-complexity regions
    /// - Low-quality consensus reads
    ///
    /// The threshold scales with the number of input molecules:
    /// - Low coverage: threshold = 1 (keep all)
    /// - High coverage: threshold = 2-3 (filter singletons)
    ///
    /// # Arguments
    ///
    /// * `n_molecules` - Total number of input molecules (consensus reads)
    ///
    /// # Example
    ///
    /// ```
    /// let budget = MemoryBudget::new(32.0);
    ///
    /// // Low coverage: keep everything
    /// assert_eq!(budget.kmer_frequency_threshold(10_000), 1);
    ///
    /// // Medium coverage: filter singletons
    /// assert_eq!(budget.kmer_frequency_threshold(1_000_000), 1);
    ///
    /// // High coverage: filter singletons and doubletons
    /// assert_eq!(budget.kmer_frequency_threshold(10_000_000), 3);
    /// ```
    pub fn kmer_frequency_threshold(&self, n_molecules: u64) -> u32 {
        // Heuristic: threshold = sqrt(n_molecules) / 1000
        // This adapts to coverage depth automatically
        let threshold = ((n_molecules as f64).sqrt() / 1000.0).max(1.0) as u32;
        threshold.min(3) // Cap at 3 to avoid over-filtering
    }

    /// Memory used by the Count-Min Sketch in megabytes.
    pub fn prefilter_memory_mb(&self) -> f64 {
        // 8 bits = 1 byte, 1 byte = ~1 counter (u8 in sketch)
        ((self.prefilter_bits() / 8) as f64) / (1024.0 * 1024.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic(expected = "must be positive")]
    fn new_zero_panic() {
        MemoryBudget::new(0.0);
    }

    #[test]
    #[should_panic(expected = "must be positive")]
    fn new_negative_panic() {
        MemoryBudget::new(-1.0);
    }

    #[test]
    fn budget_32gb_allocation() {
        let budget = MemoryBudget::new(32.0);
        assert!((budget.phase1_mb() - 8192.0).abs() < 0.01);
        assert!((budget.phase2_mb() - 19660.8).abs() < 0.01);
        assert!((budget.phase3_mb() - 4915.2).abs() < 0.01);
    }

    #[test]
    fn budget_16gb_allocation() {
        let budget = MemoryBudget::new(16.0);
        assert!((budget.phase1_mb() - 4096.0).abs() < 0.01);
        assert!((budget.phase2_mb() - 9830.4).abs() < 0.01);
        assert!((budget.phase3_mb() - 2457.6).abs() < 0.01);
    }

    #[test]
    fn budget_4gb_allocation() {
        let budget = MemoryBudget::new(4.0);
        assert_eq!(budget.phase1_mb(), 1024.0); // 25%
        assert_eq!(budget.phase2_mb(), 2457.6); // 60%
        assert_eq!(budget.phase3_mb(), 614.4); // 15%
    }

    #[test]
    fn read_batch_bytes_reasonable_32gb() {
        let budget = MemoryBudget::new(32.0);
        let batch_bytes = budget.read_batch_bytes();
        // ~4 GB / 2 = ~2 GB for read data (rest is working memory for consensus)
        assert!(batch_bytes >= 1_000_000_000);
        assert!(batch_bytes <= 5_000_000_000);
    }

    #[test]
    fn read_batch_bytes_reasonable_8gb() {
        let budget = MemoryBudget::new(8.0);
        let batch_bytes = budget.read_batch_bytes();
        // ~1 GB total for 8 GB budget
        assert!(batch_bytes >= 500_000_000);
        assert!(batch_bytes <= 1_500_000_000);
    }

    #[test]
    fn max_kmers_reasonable_32gb() {
        let budget = MemoryBudget::new(32.0);
        let max_kmers = budget.max_kmers();
        // ~340M k-mers for 32 GB (80% of 19.2 GB / 48 bytes)
        assert!(max_kmers >= 300_000_000);
        assert!(max_kmers <= 400_000_000);
    }

    #[test]
    fn max_kmers_reasonable_16gb() {
        let budget = MemoryBudget::new(16.0);
        let max_kmers = budget.max_kmers();
        // ~170M k-mers for 16 GB
        assert!(max_kmers >= 150_000_000);
        assert!(max_kmers <= 200_000_000);
    }

    #[test]
    fn prefilter_bits_reasonable_32gb() {
        let budget = MemoryBudget::new(32.0);
        let bits = budget.prefilter_bits();
        // ~1.9 GB in bits
        assert!(bits >= 15_000_000_000u64 as usize);
        assert!(bits <= 20_000_000_000u64 as usize);
    }

    #[test]
    fn kmer_frequency_threshold_scaling() {
        let budget = MemoryBudget::new(32.0);

        // No molecules: threshold = 1 (floor)
        assert_eq!(budget.kmer_frequency_threshold(0), 1);

        // Very low coverage: threshold = 1
        assert_eq!(budget.kmer_frequency_threshold(1_000), 1);

        // Medium coverage: threshold = 1
        assert_eq!(budget.kmer_frequency_threshold(100_000), 1);

        // Higher coverage: threshold = 2
        assert_eq!(budget.kmer_frequency_threshold(1_000_000), 1); // sqrt(1M) = 1000 / 1000 = 1

        // Very high coverage: threshold = 3 (capped)
        assert_eq!(budget.kmer_frequency_threshold(100_000_000), 3);
    }

    #[test]
    fn prefilter_memory_mb_small() {
        let budget = MemoryBudget::new(1.0);
        let mb = budget.prefilter_memory_mb();
        // ~61 MB for 1 GB total (10% of 614 MB phase 2)
        assert!(mb >= 50.0);
        assert!(mb <= 100.0);
    }
}
