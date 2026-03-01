//! The Disk Boundary — serialization, mmap, and columnar export.
//!
//! This layer handles bytes touching the disk: rkyv-serialized FM-Index
//! archives loaded via huge-page mmap, and zero-copy Arrow/Parquet sinks.
//! It knows how to persist and resurrect engine data structures but has
//! no concept of biology, chemistry, or Python.
//!
//! Import hierarchy: `engine`, `error` (never chemistry, operations, or lib).

pub mod parquet_sink;
pub mod persist;
