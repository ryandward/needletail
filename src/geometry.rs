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
