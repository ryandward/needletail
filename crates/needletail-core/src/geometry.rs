//! Pure coordinate geometry — modular arithmetic on 1D coordinate spaces.
//!
//! Zero biology, zero I/O, zero state. Every function is a branchless
//! mathematical operation on intervals and positions. Handles both linear
//! and circular topologies via `rem_euclid`, never via if/else boundary
//! checks.
//!
//! This module is the SeqChain coordinate algebra, ported to Rust:
//!   - `normalize`          ↔ `seqchain.primitives.coordinates.normalize`
//!   - `interval_envelope`  ↔ `seqchain.primitives.coordinates.interval_envelope`
//!   - `is_low_side`        ↔ `seqchain.primitives.coordinates.is_low_side`
//!   - `complement_base`    ↔ `seqchain.primitives.sequence.complement_base`
//!   - `fetch_sequence`     ↔ `seqchain.primitives.sequence.fetch_sequence`
//!   - `interval_overlap`   ↔ `seqchain.primitives.coordinates.interval_overlap`
//!   - `offset_in_feature`  ↔ `seqchain.primitives.coordinates.offset_in_feature`
//!   - `relative_position`  ↔ `seqchain.primitives.coordinates.relative_position`
//!   - `signed_distance`    ↔ `seqchain.primitives.coordinates.signed_distance`
//!   - `resolve_landmark`   ↔ `seqchain.primitives.coordinates.resolve_landmark`

use crate::models::preset::Anchor;

// ═══════════════════════════════════════════════════════════════════════════════
//  Coordinate arithmetic
// ═══════════════════════════════════════════════════════════════════════════════

/// Wrap a position into `[0, length)` via modular arithmetic.
/// Identity for linear (non-circular) spaces.
#[inline(always)]
pub fn normalize(pos: i64, length: i64, circular: bool) -> i64 {
    if circular {
        pos.rem_euclid(length)
    } else {
        pos
    }
}

/// Bounding envelope of two intervals. Returns `(min_start, max_end)`.
///
/// Pure min/max — no branches, no topology awareness. The caller
/// normalizes the result for circular spaces if needed.
#[inline(always)]
pub fn interval_envelope(a_start: i64, a_end: i64, b_start: i64, b_end: i64) -> (i64, i64) {
    (a_start.min(b_start), a_end.max(b_end))
}

/// Collapse a 2×2 orientation matrix into one boolean.
///
/// Given two orientation flags — `is_forward` (strand) and
/// `is_downstream` (direction) — returns `true` when the secondary
/// region sits at lower coordinates than the primary region.
///
/// Truth table:
/// ```text
///   fwd + downstream → true   (secondary below primary)
///   fwd + upstream   → false
///   rev + downstream → false
///   rev + upstream   → true
/// ```
#[inline(always)]
pub fn is_low_side(is_forward: bool, is_downstream: bool) -> bool {
    is_downstream == is_forward
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Byte-level sequence operations
// ═══════════════════════════════════════════════════════════════════════════════

/// Complement a single DNA base byte. Pure byte transform.
///
/// A↔T, C↔G (both cases). Non-ACGT → N.
#[inline(always)]
pub fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        _ => b'N',
    }
}

/// Fetch `length` bases from `text[chrom_start..]` at local offset `start`.
///
/// Circular spaces wrap positions via `rem_euclid` — a guide spanning the
/// origin is handled by modular indexing, not by boundary branches.
///
/// Linear spaces reject out-of-bounds regions (returns `false`, writes nothing).
///
/// When `rc` is true, bases are complemented and emitted 3′→5′.
#[inline(always)]
pub fn fetch_sequence(
    text: &[u8],
    chrom_start: usize,
    chrom_len: usize,
    start: i64,
    length: usize,
    circular: bool,
    rc: bool,
    dst: &mut Vec<u8>,
) -> bool {
    let end = start + length as i64;
    let cl = chrom_len as i64;
    if !circular && (start < 0 || end > cl) {
        return false;
    }
    if rc {
        for k in (0..length as i64).rev() {
            dst.push(complement_base(
                text[chrom_start + (start + k).rem_euclid(cl) as usize],
            ));
        }
    } else {
        for k in 0..length as i64 {
            dst.push(text[chrom_start + (start + k).rem_euclid(cl) as usize]);
        }
    }
    true
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Interval algebra for annotation
// ═══════════════════════════════════════════════════════════════════════════════

/// Overlap in bp between two intervals. Returns 0 if no overlap.
///
/// Linear mode: `min(a_end, b_end) - max(a_start, b_start)`, clamped to 0.
/// Circular mode: normalizes coordinates, handles wrap-around.
///
/// Port of: `seqchain.primitives.coordinates.interval_overlap`
#[inline(always)]
pub fn interval_overlap(
    a_s: i64,
    a_e: i64,
    b_s: i64,
    b_e: i64,
    chrom_len: Option<i64>,
) -> i64 {
    match chrom_len {
        None => {
            // Linear mode
            let overlap = a_e.min(b_e) - a_s.max(b_s);
            overlap.max(0)
        }
        Some(cl) => {
            // Circular mode: normalize and check overlap
            let a_s_n = a_s.rem_euclid(cl);
            let a_e_n = a_e.rem_euclid(cl);
            let b_s_n = b_s.rem_euclid(cl);
            let b_e_n = b_e.rem_euclid(cl);

            // Effective lengths
            let a_len = if a_e_n > a_s_n {
                a_e_n - a_s_n
            } else {
                cl - a_s_n + a_e_n
            };
            let b_len = if b_e_n > b_s_n {
                b_e_n - b_s_n
            } else {
                cl - b_s_n + b_e_n
            };

            // Compute gap between intervals on the circle
            let gap1 = (b_s_n - a_e_n).rem_euclid(cl);
            let gap2 = (a_s_n - b_e_n).rem_euclid(cl);
            let gap = gap1.min(gap2);

            // Total span = a_len + b_len + gap(s). If a_len + b_len > gap,
            // they overlap by the difference.
            let total = a_len + b_len;
            if total + gap <= cl {
                // They overlap by total - (cl - gap) ... actually simpler:
                // overlap = a_len + b_len - (distance traversed to cover both)
                let union = cl - gap.min(cl - a_len).min(cl - b_len);
                let _ = union; // Unused in this branch
            }

            // Simpler approach: unwrap both intervals to linear and compute
            let a_len_actual = a_e - a_s;
            let b_len_actual = b_e - b_s;
            if a_len_actual <= 0 || b_len_actual <= 0 {
                return 0;
            }
            // Try direct overlap
            let direct = a_e.min(b_e) - a_s.max(b_s);
            if direct > 0 {
                return direct;
            }
            // Try wrapping: shift b by chrom_len
            let wrap1 = a_e.min(b_e + cl) - a_s.max(b_s + cl);
            let wrap2 = (a_e + cl).min(b_e) - (a_s + cl).max(b_s);
            direct.max(wrap1).max(wrap2).max(0)
        }
    }
}

/// Strand-aware offset of a target interval within a feature.
///
/// Forward strand: `target_start - feature_start`
/// Reverse strand: `feature_end - target_end`
/// Handles circular chromosomes via modular arithmetic.
///
/// Port of: `seqchain.primitives.coordinates.offset_in_feature`
#[inline(always)]
pub fn offset_in_feature(
    t_s: i64,
    t_e: i64,
    f_s: i64,
    f_e: i64,
    fwd: bool,
    chrom_len: i64,
) -> i64 {
    if chrom_len > 0 {
        let ts = t_s.rem_euclid(chrom_len);
        let te = t_e.rem_euclid(chrom_len);
        let fs = f_s.rem_euclid(chrom_len);
        let fe = f_e.rem_euclid(chrom_len);
        if fwd {
            (ts - fs).rem_euclid(chrom_len)
        } else {
            let te_eff = if te == 0 { chrom_len } else { te };
            let fe_eff = if fe == 0 { chrom_len } else { fe };
            (fe_eff - te_eff).rem_euclid(chrom_len)
        }
    } else if fwd {
        t_s - f_s
    } else {
        f_e - t_e
    }
}

/// Fractional position within a feature. 0.0 = start, 1.0 = end.
///
/// For features that wrap the origin (virtual end > chrom_len),
/// coordinates are unwrapped before normalization.
///
/// Port of: `seqchain.primitives.coordinates.relative_position`
#[inline(always)]
pub fn relative_position(pos: i64, feat_start: i64, feat_end: i64, chrom_len: i64) -> f64 {
    let feat_len = feat_end - feat_start;
    if feat_len <= 0 {
        return 0.0;
    }

    // Handle wrapped features (virtual coordinates beyond chrom_len)
    let effective_pos = if chrom_len > 0 && feat_end > chrom_len && pos < feat_start {
        pos + chrom_len
    } else {
        pos
    };

    let offset = (effective_pos - feat_start) as f64;
    let frac = offset / feat_len as f64;
    frac.clamp(0.0, 1.0)
}

/// Strand-aware signed distance: `(query - landmark) * strand_sign`.
///
/// Positive = downstream of landmark in transcription direction.
/// Negative = upstream of landmark in transcription direction.
///
/// Port of: `seqchain.primitives.coordinates.signed_distance`
#[inline(always)]
pub fn signed_distance(query_pos: i64, landmark_pos: i64, fwd: bool) -> i64 {
    let raw = query_pos - landmark_pos;
    if fwd { raw } else { -raw }
}

/// Resolve a named landmark to a coordinate on a gene.
///
/// - FivePrime: gene start for forward, gene end for reverse (TSS)
/// - ThreePrime: gene end for forward, gene start for reverse (TES)
/// - Midpoint: (start + end) / 2 (strand-invariant)
///
/// Port of: `seqchain.primitives.coordinates.resolve_landmark`
#[inline(always)]
pub fn resolve_landmark(anchor: Anchor, gene_start: i64, gene_end: i64, fwd: bool) -> i64 {
    match anchor {
        Anchor::FivePrime => {
            if fwd { gene_start } else { gene_end }
        }
        Anchor::ThreePrime => {
            if fwd { gene_end } else { gene_start }
        }
        Anchor::Midpoint => (gene_start + gene_end) / 2,
        Anchor::None => gene_start,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interval_overlap_linear() {
        // Overlapping
        assert_eq!(interval_overlap(10, 30, 20, 40, None), 10);
        // No overlap
        assert_eq!(interval_overlap(10, 20, 30, 40, None), 0);
        // Contained
        assert_eq!(interval_overlap(10, 50, 20, 30, None), 10);
        // Same interval
        assert_eq!(interval_overlap(10, 30, 10, 30, None), 20);
    }

    #[test]
    fn test_offset_in_feature() {
        // Forward strand: target at position 15, feature starts at 10
        assert_eq!(offset_in_feature(15, 20, 10, 30, true, 0), 5);
        // Reverse strand: target ends at 25, feature ends at 30
        assert_eq!(offset_in_feature(20, 25, 10, 30, false, 0), 5);
    }

    #[test]
    fn test_relative_position() {
        assert!((relative_position(10, 0, 100, 0) - 0.1).abs() < 1e-10);
        assert!((relative_position(50, 0, 100, 0) - 0.5).abs() < 1e-10);
        assert!((relative_position(0, 0, 100, 0) - 0.0).abs() < 1e-10);
        assert!((relative_position(100, 0, 100, 0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_signed_distance() {
        // Forward strand: query downstream of landmark
        assert_eq!(signed_distance(110, 100, true), 10);
        // Forward strand: query upstream of landmark
        assert_eq!(signed_distance(90, 100, true), -10);
        // Reverse strand: query upstream is positive displacement
        assert_eq!(signed_distance(90, 100, false), 10);
    }

    #[test]
    fn test_resolve_landmark() {
        // Forward gene [100, 200)
        assert_eq!(resolve_landmark(Anchor::FivePrime, 100, 200, true), 100);
        assert_eq!(resolve_landmark(Anchor::ThreePrime, 100, 200, true), 200);
        assert_eq!(resolve_landmark(Anchor::Midpoint, 100, 200, true), 150);

        // Reverse gene [100, 200)
        assert_eq!(resolve_landmark(Anchor::FivePrime, 100, 200, false), 200);
        assert_eq!(resolve_landmark(Anchor::ThreePrime, 100, 200, false), 100);
    }
}
