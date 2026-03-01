//! CRISPR presets and feature configuration types.
//!
//! All biological knowledge enters through these configuration types.
//! The engine implements the math that configuration parameterizes.

/// Named landmark on a gene for anchor-based annotation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Anchor {
    /// Transcription start site: gene start for +, gene end for -.
    FivePrime,
    /// Transcription end site: gene end for +, gene start for -.
    ThreePrime,
    /// Gene center: (start + end) / 2.
    Midpoint,
    /// No landmark.
    None,
}

impl Anchor {
    pub fn from_str_anchor(s: &str) -> Self {
        match s {
            "five_prime" => Anchor::FivePrime,
            "three_prime" => Anchor::ThreePrime,
            "midpoint" => Anchor::Midpoint,
            _ => Anchor::None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Anchor::FivePrime => "five_prime",
            Anchor::ThreePrime => "three_prime",
            Anchor::Midpoint => "midpoint",
            Anchor::None => "",
        }
    }
}

/// A feature type defined by spatial relation to genes.
///
/// Port of SeqChain's `FeatureDefinition`.
#[derive(Debug, Clone)]
pub struct FeatureDefinition {
    /// Feature name (e.g., "promoter", "gene_body", "terminator").
    pub name: String,
    /// Spatial relation to genes: "overlap", "upstream", or "downstream".
    pub relation: String,
    /// Maximum extent in bp (0 for overlap type).
    pub max_distance: i64,
    /// Lower number = higher priority for overlap resolution.
    pub priority: i32,
    /// Named landmark for enrichment.
    pub anchor: Anchor,
}

/// Configuration for feature-type annotation.
///
/// Port of SeqChain's `FeatureConfig`.
#[derive(Debug, Clone)]
pub struct FeatureConfig {
    /// Organism name (e.g., "S. cerevisiae").
    pub organism: String,
    /// Feature definitions in priority order.
    pub features: Vec<FeatureDefinition>,
    /// Label for unmatched bases (e.g., "intergenic").
    pub default_feature: String,
}

impl FeatureConfig {
    /// Built-in S. cerevisiae feature configuration.
    pub fn saccer3() -> Self {
        FeatureConfig {
            organism: "S. cerevisiae".into(),
            features: vec![
                FeatureDefinition {
                    name: "gene_body".into(),
                    relation: "overlap".into(),
                    max_distance: 0,
                    priority: 1,
                    anchor: Anchor::FivePrime,
                },
                FeatureDefinition {
                    name: "promoter".into(),
                    relation: "upstream".into(),
                    max_distance: 500,
                    priority: 2,
                    anchor: Anchor::FivePrime,
                },
                FeatureDefinition {
                    name: "terminator".into(),
                    relation: "downstream".into(),
                    max_distance: 200,
                    priority: 3,
                    anchor: Anchor::ThreePrime,
                },
            ],
            default_feature: "intergenic".into(),
        }
    }
}

/// Configuration for a CRISPR nuclease system.
///
/// Port of SeqChain's `CRISPRPreset`.
#[derive(Debug, Clone)]
pub struct CRISPRPreset {
    pub name: String,
    pub pam: String,
    pub spacer_len: usize,
    pub pam_direction: String,
    pub description: String,
    pub mismatches: u8,
}

impl CRISPRPreset {
    /// SpCas9: NGG PAM downstream of 20bp spacer.
    pub fn spcas9() -> Self {
        CRISPRPreset {
            name: "SpCas9".into(),
            pam: "NGG".into(),
            spacer_len: 20,
            pam_direction: "downstream".into(),
            description: "Streptococcus pyogenes Cas9 (NGG PAM downstream of 20bp spacer)".into(),
            mismatches: 0,
        }
    }

    /// Cas12a (Cpf1): TTTN PAM upstream of 24bp spacer.
    pub fn cas12a() -> Self {
        CRISPRPreset {
            name: "Cas12a".into(),
            pam: "TTTN".into(),
            spacer_len: 24,
            pam_direction: "upstream".into(),
            description: "Cas12a/Cpf1 (TTTN PAM upstream of 24bp spacer)".into(),
            mismatches: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_anchor_roundtrip() {
        for anchor in [Anchor::FivePrime, Anchor::ThreePrime, Anchor::Midpoint, Anchor::None] {
            assert_eq!(Anchor::from_str_anchor(anchor.as_str()), anchor);
        }
    }

    #[test]
    fn test_spcas9_preset() {
        let p = CRISPRPreset::spcas9();
        assert_eq!(p.pam, "NGG");
        assert_eq!(p.spacer_len, 20);
        assert_eq!(p.pam_direction, "downstream");
    }

    #[test]
    fn test_saccer3_config() {
        let c = FeatureConfig::saccer3();
        assert_eq!(c.features.len(), 3);
        assert_eq!(c.features[0].name, "gene_body");
        assert_eq!(c.features[1].name, "promoter");
        assert_eq!(c.features[1].max_distance, 500);
        assert_eq!(c.default_feature, "intergenic");
    }
}
