//! Molecular Rules — the physical chemistry of CRISPR recognition.
//!
//! This module defines the payload and rules of the CRISPR enzyme:
//! PAM orientation, IUPAC bitmask compilation, and deterministic guide
//! identity hashing. It is pure computation over byte patterns — no
//! genome coordinates, no search structures, no Python.
//!
//! Import hierarchy: none (leaf module, like geometry).

use sha2::{Digest, Sha256};

// ═══════════════════════════════════════════════════════════════════════════════
//  PAM direction
// ═══════════════════════════════════════════════════════════════════════════════

/// PAM orientation relative to the spacer on the coding strand.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum PamDirection {
    Downstream,
    Upstream,
}

impl PamDirection {
    /// Convert to a raw orientation boolean for the geometry layer.
    /// `true` ≡ downstream, matching `geometry::is_low_side`'s second argument.
    #[inline(always)]
    pub fn is_downstream(self) -> bool {
        self == PamDirection::Downstream
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  IUPAC bitmask engine
// ═══════════════════════════════════════════════════════════════════════════════

/// ASCII byte → 4-bit IUPAC mask.
/// A=0b0001, C=0b0010, G=0b0100, T=0b1000. Non-ACGT → 0 (never matches).
pub static BASE_MASK: [u8; 256] = {
    let mut t = [0u8; 256];
    t[b'A' as usize] = 0b0001;
    t[b'a' as usize] = 0b0001;
    t[b'C' as usize] = 0b0010;
    t[b'c' as usize] = 0b0010;
    t[b'G' as usize] = 0b0100;
    t[b'g' as usize] = 0b0100;
    t[b'T' as usize] = 0b1000;
    t[b't' as usize] = 0b1000;
    t
};

/// Complement a 4-bit IUPAC mask: swap A↔T (bits 0↔3) and C↔G (bits 1↔2).
#[inline(always)]
fn complement_mask(m: u8) -> u8 {
    ((m & 0b0001) << 3)
        | ((m & 0b0010) << 1)
        | ((m >> 1) & 0b0010)
        | ((m >> 3) & 0b0001)
}

/// Map a single IUPAC character to its 4-bit mask.
fn iupac_to_mask(ch: u8) -> Result<u8, String> {
    match ch.to_ascii_uppercase() {
        b'A' => Ok(0b0001),
        b'C' => Ok(0b0010),
        b'G' => Ok(0b0100),
        b'T' => Ok(0b1000),
        b'N' => Ok(0b1111),
        b'R' => Ok(0b0101),
        b'Y' => Ok(0b1010),
        b'S' => Ok(0b0110),
        b'W' => Ok(0b1001),
        b'K' => Ok(0b1100),
        b'M' => Ok(0b0011),
        b'B' => Ok(0b1110),
        b'D' => Ok(0b1101),
        b'H' => Ok(0b1011),
        b'V' => Ok(0b0111),
        _ => Err(format!("invalid IUPAC code: '{}'", ch as char)),
    }
}

/// Compiled PAM pattern — forward and reverse-complement bitmask arrays.
///
/// `fwd_masks[j]` matches the PAM on the forward strand of the genome.
/// `rc_masks[j]` matches the reverse-complement PAM (complement each mask,
/// then reverse the array).
pub struct CompiledPam {
    pub fwd_masks: Vec<u8>,
    pub rc_masks: Vec<u8>,
    pub len: usize,
}

impl CompiledPam {
    pub fn compile(pam: &str) -> Result<Self, String> {
        if pam.is_empty() {
            return Err("PAM pattern must not be empty".into());
        }
        let fwd_masks: Vec<u8> = pam
            .bytes()
            .map(iupac_to_mask)
            .collect::<Result<_, _>>()?;
        let len = fwd_masks.len();
        let rc_masks: Vec<u8> = fwd_masks
            .iter()
            .rev()
            .map(|&m| complement_mask(m))
            .collect();
        Ok(CompiledPam {
            fwd_masks,
            rc_masks,
            len,
        })
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Guide identity
// ═══════════════════════════════════════════════════════════════════════════════

/// Deterministic coordinate-hashed guide identifier.
///
/// Format: 8-char hex prefix of SHA-256(`{chrom}:{start}:{strand}:{pam_seq}`).
/// Matches `seqchain.primitives.coordinates.generate_guide_id` exactly.
pub fn generate_guide_id(chrom: &str, start: u32, strand: &str, pam_seq: &str) -> [u8; 8] {
    static HEX: &[u8; 16] = b"0123456789abcdef";
    let key = format!("{}:{}:{}:{}", chrom, start, strand, pam_seq);
    let hash = Sha256::digest(key.as_bytes());
    let mut out = [0u8; 8];
    for (i, &byte) in hash[..4].iter().enumerate() {
        out[i * 2] = HEX[(byte >> 4) as usize];
        out[i * 2 + 1] = HEX[(byte & 0x0f) as usize];
    }
    out
}
