//! FmIndexSearcher — Interleaved Rank Array over a rust-bio BWT.
//!
//! Builds a BWT + suffix array from a FASTA file (via rust-bio), then
//! constructs a cache-optimized **interleaved rank array** for the hot-path
//! LF-mapping. The traditional `Occ` structure (Vec<Vec<usize>>, one heap
//! allocation per alphabet symbol) is replaced entirely.
//!
//! ## Memory layout: InterleavedRank
//!
//! ```text
//! data[pos * 4 + 0] = count of 'A' in BWT[0..=pos]
//! data[pos * 4 + 1] = count of 'C' in BWT[0..=pos]
//! data[pos * 4 + 2] = count of 'G' in BWT[0..=pos]
//! data[pos * 4 + 3] = count of 'T' in BWT[0..=pos]
//! ```
//!
//! One LF-mapping (positions `l-1` and `r`, all 4 bases) touches exactly
//! **2 cache lines** (16 bytes each, within a 64-byte line) instead of 8
//! scattered heap dereferences through bio's `Occ` structure.

use bio::alphabets::Alphabet;
use bio::data_structures::bwt::{bwt, less};
use bio::data_structures::suffix_array::{suffix_array, RawSuffixArray};
use bio::io::fasta;

use crate::error::SearchError;
use crate::simd_search::{ChromGeometry, FmOcc, BASES};

// ═══════════════════════════════════════════════════════════════════════════════
//  Interleaved Rank Array
// ═══════════════════════════════════════════════════════════════════════════════

/// Cache-optimized rank array: occ counts for A, C, G, T stored interleaved.
///
/// `data[pos * 4 + base_idx]` = count of base in `BWT[0..=pos]`.
/// A full LF-mapping (8 occ lookups: 4 bases × 2 positions) accesses exactly
/// 2 cache lines instead of 8 scattered lines in bio's per-symbol layout.
///
/// Memory: `4 * bwt_len * 4 bytes`. For SacCer3 (12M bp): ~192 MB.
/// Compare to bio's `Occ` at k=1: `85 * bwt_len * 8 bytes` ≈ 8 GB.
pub struct InterleavedRank {
    /// Flat buffer: `data[pos * 4 + base_idx]`.
    /// base_idx: A=0, C=1, G=2, T=3.
    data: Vec<u32>,
    /// Less table: `less[c]` = count of characters < c in the BWT.
    /// Indexed by raw byte value (0–255). Copied from rust-bio.
    less: [usize; 256],
}

impl InterleavedRank {
    /// Build from a BWT byte sequence and the `less` table produced by rust-bio.
    fn from_bwt_and_less(bwt_seq: &[u8], less_tbl: &[usize]) -> Self {
        let n = bwt_seq.len();
        let mut data = vec![0u32; n * 4];
        let mut counts = [0u32; 4];

        for i in 0..n {
            match bwt_seq[i] {
                b'A' => counts[0] += 1,
                b'C' => counts[1] += 1,
                b'G' => counts[2] += 1,
                b'T' => counts[3] += 1,
                _ => {} // $, N — accounted for in the less table
            }
            let base = i * 4;
            data[base]     = counts[0];
            data[base + 1] = counts[1];
            data[base + 2] = counts[2];
            data[base + 3] = counts[3];
        }

        let mut less_arr = [0usize; 256];
        for (i, &v) in less_tbl.iter().enumerate() {
            less_arr[i] = v;
        }

        InterleavedRank { data, less: less_arr }
    }

    /// Single occ lookup: count of `base` in BWT[0..=pos].
    #[inline]
    pub fn occ(&self, pos: usize, base: u8) -> u32 {
        let bi = base_idx(base);
        self.data[pos * 4 + bi]
    }

    /// Batch LF-mapping: compute all 4 `(nl, nr_exclusive)` for A, C, G, T.
    /// Loads exactly 2 contiguous 16-byte blocks — 2 cache lines.
    #[inline]
    pub fn lf_map(&self, l: usize, r: usize) -> [(u32, u32); 4] {
        let occ_l = if l > 0 {
            let b = (l - 1) * 4;
            [self.data[b], self.data[b + 1], self.data[b + 2], self.data[b + 3]]
        } else {
            [0u32; 4]
        };

        let b = r * 4;
        let occ_r = [self.data[b], self.data[b + 1], self.data[b + 2], self.data[b + 3]];

        let mut result = [(0u32, 0u32); 4];
        for bi in 0..4 {
            let less_b = self.less[BASES[bi] as usize] as u32;
            result[bi] = (less_b + occ_l[bi], less_b + occ_r[bi]);
        }
        result
    }

    /// Prefetch the cache lines needed for a future `lf_map(l, r)` call.
    /// Typically issued one (l,r) group ahead so the data arrives while the
    /// CPU processes the current group's survivor scan.
    #[inline]
    pub fn prefetch_lf(&self, l: usize, r: usize) {
        #[cfg(target_arch = "x86_64")]
        unsafe {
            use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
            if l > 0 {
                let ptr_l = self.data.as_ptr().add((l - 1) * 4) as *const i8;
                _mm_prefetch(ptr_l, _MM_HINT_T0);
            }
            let ptr_r = self.data.as_ptr().add(r * 4) as *const i8;
            _mm_prefetch(ptr_r, _MM_HINT_T0);
        }
    }

    /// Raw data slice for SIMD access.
    #[inline]
    pub fn data_slice(&self) -> &[u32] {
        &self.data
    }

    /// Less table reference for SIMD LF mapping.
    #[inline]
    pub fn less_table(&self) -> &[usize; 256] {
        &self.less
    }
}

/// Map base byte to interleaved index: A=0, C=1, G=2, T=3.
#[inline]
fn base_idx(b: u8) -> usize {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0, // unreachable in hot path; $, N never queried
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Chromosome metadata
// ═══════════════════════════════════════════════════════════════════════════════

/// Metadata for a single chromosome / FASTA record.
struct ChromInfo {
    name: String,
    /// Start offset of this chromosome in the concatenated text.
    start: usize,
    /// Length of this chromosome's sequence (excluding the '$' sentinel).
    len: usize,
}

// ═══════════════════════════════════════════════════════════════════════════════
//  FmIndexSearcher
// ═══════════════════════════════════════════════════════════════════════════════

/// FM-Index over a concatenated multi-FASTA genome.
///
/// Replaces bio's `FMIndex<BWT, Less, Occ>` with `InterleavedRank` for the
/// hot-path LF-mapping. The suffix array is kept for terminal SA lookups.
pub struct FmIndexSearcher {
    rank: InterleavedRank,
    sa: RawSuffixArray,
    chroms: Vec<ChromInfo>,
    /// Concatenated genome text (with '$' separators), kept for verify phase.
    text: Vec<u8>,
}

impl FmIndexSearcher {
    /// Build an FM-Index from a FASTA file.
    pub fn from_fasta(path: &str) -> Result<Self, SearchError> {
        let reader = fasta::Reader::from_file(path)
            .map_err(SearchError::Other)?;

        let mut text: Vec<u8> = Vec::new();
        let mut chroms: Vec<ChromInfo> = Vec::new();

        for result in reader.records() {
            let record = result.map_err(SearchError::Io)?;

            let start = text.len();
            let seq = record.seq();

            // Uppercase all bases; leave non-ACGTN as-is.
            let seq_up: Vec<u8> = seq.iter().map(|&b| b.to_ascii_uppercase()).collect();

            chroms.push(ChromInfo {
                name: record.id().to_string(),
                start,
                len: seq_up.len(),
            });

            text.extend_from_slice(&seq_up);
            // '$' (0x24) is the separator; lex-smaller than all DNA bases.
            text.push(b'$');
        }

        if chroms.is_empty() {
            return Err(SearchError::Other(anyhow::anyhow!(
                "FASTA file contains no records: {}",
                path
            )));
        }

        // Build suffix array → BWT → less table → interleaved rank array.
        //
        // The alphabet MUST include '$' so that the `less` table correctly
        // accounts for separator characters. '$' gets rank 0 (lex-smallest).
        let alphabet = Alphabet::new(b"$ACGTN");
        let sa = suffix_array(&text);
        let bwt_seq = bwt(&text, &sa);
        let less_tbl = less(&bwt_seq, &alphabet);

        // Build the interleaved rank array from the raw BWT + less table.
        // This replaces bio's Occ (Vec<Vec<usize>>) with a flat, cache-aligned
        // layout. Memory: 4 * bwt_len * 4 bytes vs 85 * bwt_len * 8 bytes.
        let rank = InterleavedRank::from_bwt_and_less(&bwt_seq, &less_tbl);

        Ok(FmIndexSearcher { rank, sa, chroms, text })
    }

    /// Chromosome names in FASTA order (index = chrom_id).
    pub fn chrom_names(&self) -> Vec<&str> {
        self.chroms.iter().map(|c| c.name.as_str()).collect()
    }

    /// Concatenated genome text for the verify phase.
    pub fn text(&self) -> &[u8] {
        &self.text
    }

    /// Extract chromosome geometry for the width-first search engine.
    pub fn chrom_geometry(&self) -> ChromGeometry {
        ChromGeometry {
            ranges: self.chroms.iter().map(|c| (c.start, c.len)).collect(),
        }
    }
}

// ─── FmOcc bridge ────────────────────────────────────────────────────────────

impl FmOcc for FmIndexSearcher {
    #[inline]
    fn less(&self, c: u8) -> usize {
        self.rank.less[c as usize]
    }

    #[inline]
    fn occ(&self, pos: usize, c: u8) -> usize {
        self.rank.occ(pos, c) as usize
    }

    #[inline]
    fn sa(&self, idx: usize) -> usize {
        self.sa[idx]
    }

    #[inline]
    fn sa_len(&self) -> usize {
        self.sa.len()
    }

    #[inline]
    fn lf_map(&self, l: usize, r: usize) -> [(u32, u32); 4] {
        self.rank.lf_map(l, r)
    }

    #[inline]
    fn prefetch_lf(&self, l: usize, r: usize) {
        self.rank.prefetch_lf(l, r);
    }

    #[inline]
    fn rank_data(&self) -> Option<(&[u32], &[usize; 256])> {
        Some((self.rank.data_slice(), self.rank.less_table()))
    }
}
