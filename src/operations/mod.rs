//! The Composition — hot loops that combine geometry, chemistry, and engine.
//!
//! This layer imports from all three foundation pillars to do the actual
//! data processing: PAM scanning, spacer extraction, off-target validation,
//! and guide enrichment. It produces flat SoA output arrays but has no
//! concept of Python or FFI.
//!
//! Import hierarchy: `geometry`, `chemistry`, `engine` (never io or lib).

pub mod pam_scanner;
