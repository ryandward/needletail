//! Sweep-line annotation — join guides against feature tiles using a
//! sorted-merge algorithm with an active window deque.
//!
//! Port of: `seqchain.operations.annotate.sweep.sorted_overlap_annotate`

use std::collections::VecDeque;

use crate::models::region::Region;

use super::locus::annotate_locus_from_features;

/// Annotate regions by sweeping them against sorted feature tiles.
///
/// Both `regions` and `features` must be sorted by (chrom, start).
/// Uses a deque-based active window: O(active_window) memory, not O(total).
///
/// For each region, finds overlapping features on the same chromosome
/// and delegates coordinate math to `annotate_locus_from_features`.
pub fn sorted_overlap_annotate(
    regions: &[Region],
    features: &[Region],
    chrom_len_fn: impl Fn(&str) -> Option<i64>,
) -> Vec<Region> {
    let mut results = Vec::with_capacity(regions.len());
    let mut active: VecDeque<usize> = VecDeque::new(); // indices into features
    let mut feat_cursor: usize = 0; // next feature to consider

    for region in regions {
        // Evict: remove features that are entirely behind the current region
        while let Some(&front) = active.front() {
            let f = &features[front];
            if f.chrom != region.chrom || f.end <= region.start {
                active.pop_front();
            } else {
                break;
            }
        }

        // Advance: load features whose start < region.end
        while feat_cursor < features.len() {
            let f = &features[feat_cursor];
            // If feature is on a later chromosome, stop
            if f.chrom > region.chrom {
                break;
            }
            // If feature starts at or beyond region end (same chrom), stop
            if f.chrom == region.chrom && f.start >= region.end {
                break;
            }
            // If feature is on the same chrom and overlaps (end > region.start),
            // add to active window
            if f.chrom == region.chrom && f.end > region.start {
                active.push_back(feat_cursor);
            }
            feat_cursor += 1;
        }

        // Filter: collect features overlapping with region on same chromosome
        let overlapping: Vec<&Region> = active
            .iter()
            .filter_map(|&idx| {
                let f = &features[idx];
                if f.chrom == region.chrom && f.start < region.end && f.end > region.start {
                    Some(f)
                } else {
                    None
                }
            })
            .collect();

        let chrom_len = chrom_len_fn(&region.chrom);
        let annotated = annotate_locus_from_features(region, &overlapping, chrom_len);
        results.extend(annotated);
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::region::Strand;

    #[test]
    fn test_sweep_annotate() {
        let guides = vec![
            Region::new("chr1", 50, 70).with_strand(Strand::Forward),
            Region::new("chr1", 150, 170).with_strand(Strand::Forward),
            Region::new("chr1", 550, 570).with_strand(Strand::Forward),
        ];

        let features = vec![
            Region::new("chr1", 0, 100)
                .with_strand(Strand::Forward)
                .with_name("promoter_A")
                .with_tag("feature_type", "promoter"),
            Region::new("chr1", 100, 500)
                .with_strand(Strand::Forward)
                .with_name("gene_A")
                .with_tag("feature_type", "gene_body"),
            Region::new("chr1", 500, 600)
                .with_strand(Strand::Forward)
                .with_name("term_A")
                .with_tag("feature_type", "terminator"),
        ];

        let results = sorted_overlap_annotate(&guides, &features, |_| None);

        // Each guide should get annotated against one feature
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].tags["feature_name"].as_str(), Some("promoter_A"));
        assert_eq!(results[1].tags["feature_name"].as_str(), Some("gene_A"));
        assert_eq!(results[2].tags["feature_name"].as_str(), Some("term_A"));
    }

    #[test]
    fn test_sweep_no_overlap() {
        let guides = vec![
            Region::new("chr1", 1000, 1020).with_strand(Strand::Forward),
        ];
        let features = vec![
            Region::new("chr1", 0, 100)
                .with_name("feat")
                .with_tag("feature_type", "gene"),
        ];

        let results = sorted_overlap_annotate(&guides, &features, |_| None);
        assert_eq!(results.len(), 1);
        // No feature_name tag — guide returned unchanged
        assert!(!results[0].tags.contains_key("feature_name"));
    }
}
