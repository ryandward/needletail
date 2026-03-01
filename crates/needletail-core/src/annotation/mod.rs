//! Layer 2 — Interval algorithms for genome annotation.
//!
//! Depends on `geometry.rs` for coordinate math and `models/` for Region.
//! Parallel to `engine/` — shares geometry but doesn't know about search.

pub mod feature;
pub mod locus;
pub mod sweep;
