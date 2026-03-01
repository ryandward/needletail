//! Layer 0 — Pure data types with no logic and no dependencies on other modules.
//!
//! Region is the universal currency: every pipeline stage consumes and produces
//! Regions. Tags carry domain-specific metadata without encoding it in the type
//! system.

pub mod genome;
pub mod preset;
pub mod region;

pub use genome::{Genome, Topology};
pub use preset::{Anchor, CRISPRPreset, FeatureConfig, FeatureDefinition};
pub use region::{Region, Strand, TagValue};
