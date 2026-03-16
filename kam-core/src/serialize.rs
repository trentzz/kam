//! Bincode serialization for inter-stage data handoff.
//!
//! Every kam bincode file starts with a [`FileHeader`] that identifies the
//! file type and format version, followed by a sequence of records encoded
//! with bincode.  This allows `kam explore` to inspect any intermediate file
//! without knowing which stage produced it.
//!
//! # Example
//! ```no_run
//! use std::path::Path;
//! use kam_core::serialize::{FileType, write_bincode, read_bincode};
//!
//! let records = vec![1_u32, 2, 3];
//! let path = Path::new("/tmp/example.bin");
//! write_bincode(path, FileType::Molecules, &records).unwrap();
//! let (header, loaded): (_, Vec<u32>) = read_bincode(path).unwrap();
//! assert_eq!(loaded, records);
//! ```

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use serde::{Deserialize, Serialize};

/// File type tag written in the header of every kam bincode file.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum FileType {
    /// A sequence of [`crate::molecule::Molecule`] records.
    Molecules,
    /// A k-mer index (map of `u64` → [`crate::kmer::MoleculeEvidence`]).
    KmerIndex,
    /// Scored de Bruijn paths.
    ScoredPaths,
    /// Variant call records.
    VariantCalls,
}

/// Header written at the start of every kam bincode file.
///
/// The magic bytes `b"KAM\0"` let downstream tools quickly verify they are
/// reading a kam-produced file.  [`file_type`](FileHeader::file_type) and
/// [`version`](FileHeader::version) guard against type confusion and format
/// drift.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileHeader {
    /// Magic bytes — always `b"KAM\0"`.
    pub magic: [u8; 4],
    /// Format version (currently `1`).
    pub version: u32,
    /// Which kind of records follow the header.
    pub file_type: FileType,
    /// Number of records written after the header.
    pub record_count: u64,
}

impl FileHeader {
    /// The magic bytes used in every kam bincode file.
    pub const MAGIC: [u8; 4] = *b"KAM\0";
    /// Current format version.
    pub const VERSION: u32 = 1;

    fn new(file_type: FileType, record_count: u64) -> Self {
        Self {
            magic: Self::MAGIC,
            version: Self::VERSION,
            file_type,
            record_count,
        }
    }
}

/// Write a header followed by `records` to a bincode file at `path`.
///
/// The file is created (or truncated) and wrapped in a [`BufWriter`] for
/// efficiency.  The header [`record_count`](FileHeader::record_count) is set
/// to `records.len()`.
///
/// # Errors
///
/// Returns an error if the file cannot be created or if bincode serialization
/// fails.
///
/// # Example
/// ```no_run
/// use std::path::Path;
/// use kam_core::serialize::{FileType, write_bincode};
///
/// write_bincode(Path::new("/tmp/out.bin"), FileType::Molecules, &[42_u64]).unwrap();
/// ```
pub fn write_bincode<T: Serialize>(
    path: &Path,
    file_type: FileType,
    records: &[T],
) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    let header = FileHeader::new(file_type, records.len() as u64);
    let header_bytes = bincode::serialize(&header)?;
    writer.write_all(&header_bytes)?;

    for record in records {
        let record_bytes = bincode::serialize(record)?;
        writer.write_all(&record_bytes)?;
    }

    writer.flush()?;
    Ok(())
}

/// Read only the header from a bincode file, without loading any records.
///
/// This is useful for quickly inspecting a file's type and record count
/// without the memory cost of deserializing all records.
///
/// # Errors
///
/// Returns an error if the file cannot be opened, if the magic bytes are
/// wrong, or if the header cannot be deserialized.
///
/// # Example
/// ```no_run
/// use std::path::Path;
/// use kam_core::serialize::{FileType, write_bincode, read_header};
///
/// write_bincode(Path::new("/tmp/out.bin"), FileType::KmerIndex, &[1_u32, 2, 3]).unwrap();
/// let hdr = read_header(Path::new("/tmp/out.bin")).unwrap();
/// assert_eq!(hdr.record_count, 3);
/// ```
pub fn read_header(path: &Path) -> Result<FileHeader, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);

    // Read the entire file into memory so bincode can decode the header.
    // We only need the first few bytes, but bincode requires a contiguous
    // buffer — reading a fixed-size prefix is the simplest safe approach.
    let mut buf = Vec::new();
    reader.read_to_end(&mut buf)?;

    let header: FileHeader = bincode::deserialize_from(&mut buf.as_slice())?;

    if header.magic != FileHeader::MAGIC {
        return Err(format!(
            "invalid magic bytes: expected {:?}, got {:?}",
            FileHeader::MAGIC,
            header.magic
        )
        .into());
    }

    Ok(header)
}

/// Read the header and all records from a bincode file.
///
/// Returns `(header, records)`.  The type parameter `T` must match the record
/// type that was written; using the wrong type will produce a deserialization
/// error.
///
/// # Errors
///
/// Returns an error if the file cannot be opened, if the magic bytes are
/// wrong, or if any record cannot be deserialized.
///
/// # Example
/// ```no_run
/// use std::path::Path;
/// use kam_core::serialize::{FileType, write_bincode, read_bincode};
///
/// let data = vec![10_u32, 20, 30];
/// write_bincode(Path::new("/tmp/out.bin"), FileType::Molecules, &data).unwrap();
/// let (hdr, loaded): (_, Vec<u32>) = read_bincode(Path::new("/tmp/out.bin")).unwrap();
/// assert_eq!(loaded, data);
/// ```
pub fn read_bincode<T: for<'de> Deserialize<'de>>(
    path: &Path,
) -> Result<(FileHeader, Vec<T>), Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);

    let mut buf = Vec::new();
    reader.read_to_end(&mut buf)?;

    let mut cursor = std::io::Cursor::new(&buf);

    let header: FileHeader = bincode::deserialize_from(&mut cursor)?;

    if header.magic != FileHeader::MAGIC {
        return Err(format!(
            "invalid magic bytes: expected {:?}, got {:?}",
            FileHeader::MAGIC,
            header.magic
        )
        .into());
    }

    let mut records = Vec::with_capacity(header.record_count as usize);
    for _ in 0..header.record_count {
        let record: T = bincode::deserialize_from(&mut cursor)?;
        records.push(record);
    }

    Ok((header, records))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Seek, SeekFrom};

    fn tmp_path(name: &str) -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("kam_core_serialize_test_{name}_{}", std::process::id()));
        p
    }

    // 1. Write + read roundtrip produces identical data.
    #[test]
    fn roundtrip_produces_identical_data() {
        let path = tmp_path("roundtrip");
        let records = vec![10_u32, 20, 30, 40];
        write_bincode(&path, FileType::Molecules, &records).unwrap();
        let (_, loaded): (_, Vec<u32>) = read_bincode(&path).unwrap();
        assert_eq!(loaded, records);
        let _ = std::fs::remove_file(&path);
    }

    // 2. Header magic bytes are correct.
    #[test]
    fn header_magic_bytes_correct() {
        let path = tmp_path("magic");
        write_bincode(&path, FileType::KmerIndex, &[1_u8]).unwrap();
        let header = read_header(&path).unwrap();
        assert_eq!(header.magic, *b"KAM\0");
        let _ = std::fs::remove_file(&path);
    }

    // 3. Header record_count matches input.
    #[test]
    fn header_record_count_matches_input() {
        let path = tmp_path("count");
        let records: Vec<u64> = vec![1, 2, 3, 4, 5];
        write_bincode(&path, FileType::ScoredPaths, &records).unwrap();
        let header = read_header(&path).unwrap();
        assert_eq!(header.record_count, 5);
        let _ = std::fs::remove_file(&path);
    }

    // 4. read_header alone works without loading records.
    #[test]
    fn read_header_alone_without_loading_records() {
        let path = tmp_path("header_only");
        let records = vec![42_u32; 100];
        write_bincode(&path, FileType::VariantCalls, &records).unwrap();
        // Should succeed and return the correct type without deserialising records.
        let header = read_header(&path).unwrap();
        assert_eq!(header.file_type, FileType::VariantCalls);
        assert_eq!(header.record_count, 100);
        let _ = std::fs::remove_file(&path);
    }

    // 5. Wrong magic bytes returns error.
    #[test]
    fn wrong_magic_bytes_returns_error() {
        let path = tmp_path("bad_magic");

        // Write a valid file first, then corrupt the first 4 bytes.
        write_bincode(&path, FileType::Molecules, &[0_u8]).unwrap();
        {
            // Open for writing at offset 0 to overwrite the magic bytes.
            let mut f = std::fs::OpenOptions::new().write(true).open(&path).unwrap();
            f.seek(SeekFrom::Start(0)).unwrap();
            f.write_all(b"BAD!").unwrap();
        }

        let result = read_header(&path);
        assert!(result.is_err(), "expected error for bad magic bytes");
        let _ = std::fs::remove_file(&path);
    }

    // 6. Empty records list writes a valid file with record_count 0.
    #[test]
    fn empty_records_writes_valid_file_with_count_zero() {
        let path = tmp_path("empty");
        let records: Vec<u32> = vec![];
        write_bincode(&path, FileType::Molecules, &records).unwrap();
        let header = read_header(&path).unwrap();
        assert_eq!(header.record_count, 0);
        let (_, loaded): (_, Vec<u32>) = read_bincode(&path).unwrap();
        assert!(loaded.is_empty());
        let _ = std::fs::remove_file(&path);
    }
}
