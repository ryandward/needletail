//! Genome — named collection of sequences, features, and topologies.

use std::collections::HashMap;

use super::region::Region;

/// Chromosome topology.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Topology {
    Linear,
    Circular,
}

impl Topology {
    pub fn is_circular(self) -> bool {
        self == Topology::Circular
    }

    pub fn from_str_topo(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "circular" => Topology::Circular,
            _ => Topology::Linear,
        }
    }
}

/// A genome: sequences, features, and metadata for one or more chromosomes.
#[derive(Debug, Clone)]
pub struct Genome {
    pub name: String,
    /// Chromosome name → DNA sequence (uppercase bytes). Insertion-ordered.
    pub sequences: Vec<(String, Vec<u8>)>,
    /// All features across all chromosomes.
    pub features: Vec<Region>,
    /// Chromosome name → topology.
    pub topologies: HashMap<String, Topology>,
}

impl Genome {
    pub fn new(name: impl Into<String>) -> Self {
        Genome {
            name: name.into(),
            sequences: Vec::new(),
            features: Vec::new(),
            topologies: HashMap::new(),
        }
    }

    /// Chromosome lengths in insertion order.
    pub fn chrom_lengths(&self) -> Vec<(&str, usize)> {
        self.sequences
            .iter()
            .map(|(name, seq)| (name.as_str(), seq.len()))
            .collect()
    }

    /// Chromosome lengths as a HashMap.
    pub fn chrom_length_map(&self) -> HashMap<&str, usize> {
        self.sequences
            .iter()
            .map(|(name, seq)| (name.as_str(), seq.len()))
            .collect()
    }

    /// Sorted list of chromosome names.
    pub fn chroms(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.sequences.iter().map(|(n, _)| n.as_str()).collect();
        names.sort();
        names
    }

    /// Filter features to genes only (tags["feature_type"] == "gene").
    pub fn genes(&self) -> Vec<&Region> {
        self.features
            .iter()
            .filter(|r| {
                r.tags
                    .get("feature_type")
                    .and_then(|v| v.as_str())
                    .map_or(false, |t| t == "gene")
            })
            .collect()
    }

    /// Filter features to a single chromosome.
    pub fn features_on(&self, chrom: &str) -> Vec<&Region> {
        self.features.iter().filter(|r| r.chrom == chrom).collect()
    }

    /// Check if a chromosome has circular topology.
    pub fn is_circular(&self, chrom: &str) -> bool {
        self.topologies
            .get(chrom)
            .map_or(false, |t| t.is_circular())
    }

    /// Get sequence for a chromosome.
    pub fn sequence(&self, chrom: &str) -> Option<&[u8]> {
        self.sequences
            .iter()
            .find(|(n, _)| n == chrom)
            .map(|(_, s)| s.as_slice())
    }

    /// Boolean topology vector in insertion order (for compatibility with engine).
    pub fn topology_vec(&self) -> Vec<bool> {
        self.sequences
            .iter()
            .map(|(name, _)| self.is_circular(name))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::region::{Region, Strand};

    #[test]
    fn test_genome_basics() {
        let mut g = Genome::new("test");
        g.sequences.push(("chr1".into(), b"ATCGATCG".to_vec()));
        g.sequences.push(("chr2".into(), b"GGCC".to_vec()));
        g.topologies.insert("chr1".into(), Topology::Linear);
        g.topologies.insert("chr2".into(), Topology::Circular);

        assert_eq!(g.chrom_lengths(), vec![("chr1", 8), ("chr2", 4)]);
        assert!(!g.is_circular("chr1"));
        assert!(g.is_circular("chr2"));
        assert_eq!(g.sequence("chr1"), Some(b"ATCGATCG".as_slice()));
    }

    #[test]
    fn test_genome_genes() {
        let mut g = Genome::new("test");
        g.sequences.push(("chr1".into(), b"ATCGATCG".to_vec()));
        g.features.push(
            Region::new("chr1", 0, 100)
                .with_strand(Strand::Forward)
                .with_tag("feature_type", "gene"),
        );
        g.features.push(
            Region::new("chr1", 200, 300)
                .with_strand(Strand::Forward)
                .with_tag("feature_type", "CDS"),
        );

        let genes = g.genes();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].start, 0);
    }
}
