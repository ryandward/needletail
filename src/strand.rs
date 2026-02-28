//! StrandProducer — Rust mirror of the @strand/core binary protocol.
//!
//! Constants mirror strand/src/constants.ts verbatim; any change there is a
//! breaking change requiring a version bump in STRAND_VERSION.
//!
//! Safety invariant: The raw pointer stored in `index_base` must remain valid
//! for the entire lifetime of the StrandProducer. In practice the caller holds
//! a SharedArrayBuffer that owns the memory; it must not be GC'd during use.
//!
//! Note on AtomicU32 vs AtomicI32:
//!   The Strand protocol specifies Int32Array control words. Rust uses AtomicU32
//!   here because the `atomic_wait` crate only supports AtomicU32. The control
//!   word values (sequence counters, status, abort flag) are small positive
//!   integers whose bit patterns are identical whether interpreted as i32 or u32.
//!   JS Atomics.store/load on Int32Array and Rust AtomicU32 are futex-compatible
//!   on the same memory address — both use the same 4-byte storage cell.

use std::ptr::NonNull;
use std::sync::atomic::{AtomicU32, Ordering};
use std::time::Duration;

use crate::error::StrandError;

// ─── Header constants (mirror strand/src/constants.ts) ───────────────────────

const HEADER_SIZE: usize = 512;
const STRAND_MAGIC: u32 = 0x5354_524E; // 'STRN' little-endian
const STRAND_VERSION: u32 = 4; // v4: added producer metadata region in header tail

// Int32Array indices (i × 4 = byte offset within SAB)
#[allow(dead_code)]
const CTRL_WRITE_SEQ: usize = 7; // byte 28 — not written by Rust producer (uses local counter)
const CTRL_COMMIT_SEQ: usize = 8; // byte 32
const CTRL_READ_CURSOR: usize = 9; // byte 36
const CTRL_STATUS: usize = 10; // byte 40
const CTRL_ABORT: usize = 13; // byte 52

// Header byte offsets for static fields (DataView reads)
const OFFSET_MAGIC: usize = 0;
const OFFSET_VERSION: usize = 4;
const OFFSET_SCHEMA_FP: usize = 8;
const OFFSET_RECORD_STRIDE: usize = 12;
const OFFSET_INDEX_CAPACITY: usize = 16;
const OFFSET_HEADER_CRC: usize = 24;
const OFFSET_SCHEMA_BYTE_LEN: usize = 88;
const OFFSET_SCHEMA_BYTES: usize = 92;

// Status values (match STATUS_* in constants.ts)
const STATUS_STREAMING: u32 = 1;
const STATUS_EOS: u32 = 2;
const STATUS_ERROR: u32 = 3;

// Schema byte offsets produced by buildSchema() for our alignment schema.
//
// buildSchema([
//   { name: 'chrom',      type: 'utf8_ref' },  // width 4, align 4 → offset  0
//   { name: 'pos',        type: 'u32'      },  // width 4, align 4 → offset  4
//   { name: 'query_id',   type: 'u32'      },  // width 4, align 4 → offset  8
//   { name: 'strand',     type: 'bool8'    },  // width 1, align 1 → offset 12
//   { name: 'mismatches', type: 'u8'       },  // width 1, align 1 → offset 13
//   { name: 'score',      type: 'f32'      },  // width 4, align 4 → offset 16
// ])                                            // (padded 2 bytes from 14 to align f32)
// record_stride = 20
const OFF_CHROM: usize = 0;
const OFF_POS: usize = 4;
const OFF_QUERY_ID: usize = 8;
const OFF_STRAND: usize = 12;
const OFF_MISMATCHES: usize = 13;
const OFF_SCORE: usize = 16;
pub(crate) const RECORD_STRIDE: usize = 20;

// ─── FNV-1a 32-bit (mirrors the JS implementation in schema.ts) ──────────────

fn fnv1a32(bytes: &[u8]) -> u32 {
    let mut hash: u32 = 0x811c_9dc5;
    for &byte in bytes {
        hash ^= byte as u32;
        hash = hash.wrapping_mul(0x0100_0193);
    }
    hash
}

/// Compute the header geometry CRC: FNV-1a of bytes 0–23.
fn compute_header_crc(sab: &[u8]) -> u32 {
    fnv1a32(&sab[0..OFFSET_HEADER_CRC])
}

// ─── StrandProducer ──────────────────────────────────────────────────────────

/// Zero-copy producer that writes alignment records directly into a
/// Strand SharedArrayBuffer's index ring.
pub struct StrandProducer {
    /// Pointer to the AtomicU32 control words (128 entries = 512 bytes).
    ctrl: NonNull<[AtomicU32; 128]>,
    /// Pointer to byte 512 of the SAB (start of the index ring).
    index_base: NonNull<u8>,
    /// Number of ring slots (always a power of two).
    index_capacity: usize,
    /// Producer's local monotonic write sequence counter.
    write_seq: u32,
}

impl StrandProducer {
    /// Validate the Strand header and construct a producer.
    ///
    /// # Safety
    ///
    /// `ptr` must point to at least `len` bytes of SAB memory that remains
    /// alive and accessible for the duration of this producer's lifetime.
    pub unsafe fn new(ptr: *mut u8, len: usize) -> Result<Self, StrandError> {
        if len < HEADER_SIZE {
            return Err(StrandError::InvalidHeader(format!(
                "SAB too small: {} bytes, need at least {}",
                len, HEADER_SIZE
            )));
        }

        let sab = std::slice::from_raw_parts(ptr, len);

        // ── Magic ────────────────────────────────────────────────────────────
        let magic = read_u32_le(sab, OFFSET_MAGIC);
        if magic != STRAND_MAGIC {
            return Err(StrandError::InvalidHeader(format!(
                "Bad magic: 0x{:08X}, expected 0x{:08X}",
                magic, STRAND_MAGIC
            )));
        }

        // ── Version ──────────────────────────────────────────────────────────
        let version = read_u32_le(sab, OFFSET_VERSION);
        if version != STRAND_VERSION {
            return Err(StrandError::InvalidHeader(format!(
                "Unsupported version {}, expected {}",
                version, STRAND_VERSION
            )));
        }

        // ── Header CRC (FNV-1a of bytes 0–23) ────────────────────────────────
        let stored_crc = read_u32_le(sab, OFFSET_HEADER_CRC);
        let computed_crc = compute_header_crc(sab);
        if stored_crc != computed_crc {
            return Err(StrandError::InvalidHeader(format!(
                "Header CRC mismatch: stored 0x{:08X}, computed 0x{:08X}",
                stored_crc, computed_crc
            )));
        }

        // ── Geometry ─────────────────────────────────────────────────────────
        let index_capacity = read_u32_le(sab, OFFSET_INDEX_CAPACITY) as usize;
        if index_capacity == 0 || (index_capacity & (index_capacity - 1)) != 0 {
            return Err(StrandError::InvalidHeader(format!(
                "index_capacity {} is not a power of two",
                index_capacity
            )));
        }

        let record_stride = read_u32_le(sab, OFFSET_RECORD_STRIDE) as usize;
        if record_stride != RECORD_STRIDE {
            return Err(StrandError::InvalidHeader(format!(
                "record_stride {} != expected {}",
                record_stride, RECORD_STRIDE
            )));
        }

        // ── Schema fingerprint ────────────────────────────────────────────────
        // Re-hash the schema bytes from the header and compare to stored fp.
        let stored_fp = read_u32_le(sab, OFFSET_SCHEMA_FP);
        let schema_byte_len = read_u32_le(sab, OFFSET_SCHEMA_BYTE_LEN) as usize;
        if OFFSET_SCHEMA_BYTES + schema_byte_len > len {
            return Err(StrandError::InvalidHeader(
                "schema_byte_len extends beyond SAB bounds".to_string(),
            ));
        }
        let schema_bytes = &sab[OFFSET_SCHEMA_BYTES..OFFSET_SCHEMA_BYTES + schema_byte_len];
        let computed_fp = fnv1a32(schema_bytes);
        if stored_fp != computed_fp {
            return Err(StrandError::InvalidHeader(format!(
                "schema_fp mismatch: stored 0x{:08X}, computed 0x{:08X}",
                stored_fp, computed_fp
            )));
        }

        // ── Construct ─────────────────────────────────────────────────────────
        let ctrl = NonNull::new_unchecked(ptr as *mut [AtomicU32; 128]);
        let index_base = NonNull::new_unchecked(ptr.add(HEADER_SIZE));

        Ok(StrandProducer {
            ctrl,
            index_base,
            index_capacity,
            write_seq: 0,
        })
    }

    /// Signal STREAMING status. Call before the first `write_record`.
    pub fn begin(&self) {
        self.ctrl()[CTRL_STATUS].store(STATUS_STREAMING, Ordering::Release);
    }

    /// Write one alignment record into the ring buffer.
    ///
    /// Blocks (via Linux futex / atomic-wait) if the ring is full.
    /// Returns `Err(Aborted)` immediately if the consumer signals abort.
    pub fn write_record(
        &mut self,
        chrom: u32,    // utf8_ref handle — index into the intern table
        pos: u32,
        query_id: u32,
        strand: bool,
        mismatches: u8,
        score: f32,
    ) -> Result<(), StrandError> {
        // Use raw pointer to avoid borrow-checker conflicts between ctrl reads
        // (shared borrow of self) and write_seq mutation (mut borrow of self).
        // Safety: `self.ctrl` points into the SAB which outlives this call.
        let ctrl_ptr: *const [AtomicU32; 128] = self.ctrl.as_ptr();
        let ctrl = unsafe { &*ctrl_ptr };

        // ── Abort check ───────────────────────────────────────────────────────
        if ctrl[CTRL_ABORT].load(Ordering::Acquire) != 0 {
            return Err(StrandError::Aborted);
        }

        // ── Backpressure stall ────────────────────────────────────────────────
        // Block until write_seq - read_cursor < index_capacity.
        //
        // Note: V8's Atomics.notify uses an internal user-space futex
        // (condition variables), while atomic_wait::wait uses the Linux kernel
        // futex syscall. They don't interoperate. We use a short spin-sleep
        // instead, which reacts to READ_CURSOR changes within ~50µs.
        let capacity = self.index_capacity as u32;
        loop {
            let read_cursor = ctrl[CTRL_READ_CURSOR].load(Ordering::Acquire);
            if self.write_seq.wrapping_sub(read_cursor) < capacity {
                break;
            }
            if ctrl[CTRL_ABORT].load(Ordering::Acquire) != 0 {
                return Err(StrandError::Aborted);
            }
            // Short sleep: yields the thread and re-checks every 50µs.
            // The JS consumer advances READ_CURSOR via acknowledgeRead().
            std::thread::sleep(Duration::from_micros(50));
        }

        // ── Compute slot and base pointer ─────────────────────────────────────
        let slot = (self.write_seq as usize) & (self.index_capacity - 1);
        let record_ptr = unsafe { self.index_base.as_ptr().add(slot * RECORD_STRIDE) };

        // ── Write fields ──────────────────────────────────────────────────────
        unsafe {
            write_u32_le(record_ptr.add(OFF_CHROM), chrom);
            write_u32_le(record_ptr.add(OFF_POS), pos);
            write_u32_le(record_ptr.add(OFF_QUERY_ID), query_id);
            record_ptr.add(OFF_STRAND).write(strand as u8);
            record_ptr.add(OFF_MISMATCHES).write(mismatches);
            write_f32_le(record_ptr.add(OFF_SCORE), score);
        }

        // ── Advance and commit ────────────────────────────────────────────────
        // Mutation of write_seq is safe here: ctrl is a raw-pointer-derived ref
        // that does not alias write_seq (different SAB region vs struct field).
        self.write_seq = self.write_seq.wrapping_add(1);
        ctrl[CTRL_COMMIT_SEQ].store(self.write_seq, Ordering::Release);
        // Note: atomic_wait::wake_one uses kernel futex WAKE, which does not
        // interoperate with V8's Atomics.waitAsync (user-space futex). The JS
        // consumer should use a polling loop with a timeout rather than
        // Atomics.waitAsync for zero-latency consumption.

        Ok(())
    }

    /// Signal EOS and wake the consumer. Call after all records are written.
    pub fn finalize(&self) {
        self.ctrl()[CTRL_STATUS].store(STATUS_EOS, Ordering::Release);
        // JS consumer detects EOS via status poll or Promise resolution.
    }

    /// Signal ERROR and wake the consumer. Call on panic or write error.
    pub fn abort(&self) {
        self.ctrl()[CTRL_STATUS].store(STATUS_ERROR, Ordering::Release);
    }

    /// Access the AtomicU32 control word array.
    fn ctrl(&self) -> &[AtomicU32; 128] {
        unsafe { self.ctrl.as_ref() }
    }
}

// ─── Field write helpers ──────────────────────────────────────────────────────

#[inline(always)]
unsafe fn write_u32_le(ptr: *mut u8, val: u32) {
    ptr.copy_from_nonoverlapping(val.to_le_bytes().as_ptr(), 4);
}

#[inline(always)]
unsafe fn write_f32_le(ptr: *mut u8, val: f32) {
    ptr.copy_from_nonoverlapping(val.to_le_bytes().as_ptr(), 4);
}

// ─── Header read helpers ──────────────────────────────────────────────────────

#[inline(always)]
fn read_u32_le(sab: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(sab[offset..offset + 4].try_into().unwrap())
}
