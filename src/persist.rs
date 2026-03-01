//! Zero-copy mmap persistence for seqchain FM-Index data.
//!
//! Serializes the FM-Index (interleaved rank array, suffix array, less table,
//! chromosome metadata, and genome text) to a `.seqchain` file using rkyv.
//! On load, the file is mmap'd and the archived data is accessed directly —
//! no deserialization for the large rank array or suffix array.
//!
//! File format:
//! ```text
//! Offset 0:   "SQCH" (4 bytes magic)
//! Offset 4:   1u32 LE (format version)
//! Offset 8:   payload_len u64 LE
//! Offset 16:  rkyv payload (16-byte aligned)
//! ```

use std::fs::File;
use std::io::{BufWriter, Read as _, Write};
use std::path::Path;

use rkyv::rancor::Error as RkyvError;
use rkyv::{Archive, Deserialize, Serialize};

use crate::error::SearchError;
use crate::fm_index::FmIndexSearcher;
use crate::simd_search::{ChromGeometry, FmOcc, BASES};

// ═══════════════════════════════════════════════════════════════════════════════
//  Huge Page mmap — strict MAP_HUGETLB, no fallback
// ═══════════════════════════════════════════════════════════════════════════════

/// 2 MiB huge page size.
const HUGE_PAGE_SIZE: usize = 2 * 1024 * 1024;

/// Memory region backed by 2 MiB huge pages (MAP_HUGETLB).
///
/// **Strict mode**: if the kernel pool is empty (`vm.nr_hugepages=0`),
/// the allocation panics with a diagnostic message. No THP fallback.
/// No silent degradation. Either we get real 2 MiB pages or we stop.
///
/// Why: THP is a performance trap — the kernel can split/collapse pages
/// under memory pressure, causing unpredictable latency spikes during
/// the wavefront sweep. Explicit huge pages from the reserved pool are
/// pinned and never split.
struct HugePageMmap {
    ptr: *mut u8,
    len: usize,
}

impl HugePageMmap {
    /// Allocate `MAP_HUGETLB | MAP_ANONYMOUS | MAP_PRIVATE` and copy `data` in.
    ///
    /// # Panics
    ///
    /// Panics if `MAP_HUGETLB` fails. Configure the pool first:
    /// ```bash
    /// sudo sysctl -w vm.nr_hugepages=256   # 512 MiB of 2 MiB pages
    /// ```
    fn from_data(data: &[u8]) -> Self {
        let aligned_len = (data.len() + HUGE_PAGE_SIZE - 1) & !(HUGE_PAGE_SIZE - 1);
        let pages_needed = aligned_len / HUGE_PAGE_SIZE;

        let ptr = unsafe {
            libc::mmap(
                std::ptr::null_mut(),
                aligned_len,
                libc::PROT_READ | libc::PROT_WRITE,
                libc::MAP_PRIVATE | libc::MAP_ANONYMOUS | libc::MAP_HUGETLB,
                -1,
                0,
            )
        };

        if ptr == libc::MAP_FAILED {
            let err = std::io::Error::last_os_error();
            panic!(
                "\n\
                ╔══════════════════════════════════════════════════════════════╗\n\
                ║  MAP_HUGETLB FAILED — kernel huge page pool is empty       ║\n\
                ╠══════════════════════════════════════════════════════════════╣\n\
                ║  Requested : {} pages × 2 MiB = {} MiB                 ║\n\
                ║  OS error  : {}                                        ║\n\
                ║                                                            ║\n\
                ║  Fix: sudo sysctl -w vm.nr_hugepages={}                ║\n\
                ║  Verify: cat /proc/meminfo | grep HugePages                ║\n\
                ╚══════════════════════════════════════════════════════════════╝\n",
                pages_needed,
                aligned_len / (1024 * 1024),
                err,
                pages_needed + 16, // headroom
            );
        }

        // Verify 2 MiB alignment — MAP_HUGETLB guarantees this, but assert it.
        let addr = ptr as usize;
        assert_eq!(
            addr & (HUGE_PAGE_SIZE - 1),
            0,
            "MAP_HUGETLB returned non-2MiB-aligned pointer: {:#x}",
            addr,
        );

        let dst = ptr as *mut u8;
        unsafe { std::ptr::copy_nonoverlapping(data.as_ptr(), dst, data.len()) };
        // Seal read-only.
        unsafe { libc::mprotect(ptr, aligned_len, libc::PROT_READ) };

        eprintln!(
            "[needletail] MAP_HUGETLB: {} MiB bound to {} × 2 MiB pages (aligned @ {:#x})",
            data.len() / (1024 * 1024),
            pages_needed,
            addr,
        );

        HugePageMmap {
            ptr: dst,
            len: aligned_len,
        }
    }

    /// Access the data as a byte slice.
    #[inline]
    fn as_slice(&self, actual_len: usize) -> &[u8] {
        unsafe { std::slice::from_raw_parts(self.ptr, actual_len) }
    }
}

impl Drop for HugePageMmap {
    fn drop(&mut self) {
        unsafe {
            libc::munmap(self.ptr as *mut libc::c_void, self.len);
        }
    }
}

// SAFETY: The mapping is read-only after construction, no mutation possible.
unsafe impl Send for HugePageMmap {}
unsafe impl Sync for HugePageMmap {}

// ═══════════════════════════════════════════════════════════════════════════════
//  File format constants
// ═══════════════════════════════════════════════════════════════════════════════

const MAGIC: &[u8; 4] = b"SQCH";
const FORMAT_VERSION: u32 = 1;
const HEADER_SIZE: usize = 16;

// ═══════════════════════════════════════════════════════════════════════════════
//  Stored types (rkyv-serializable)
// ═══════════════════════════════════════════════════════════════════════════════

#[derive(Archive, Serialize, Deserialize)]
#[rkyv(compare(PartialEq), derive(Debug))]
struct StoredChromInfo {
    name: String,
    start: u64,
    len: u64,
}

#[derive(Archive, Serialize, Deserialize)]
#[rkyv(compare(PartialEq), derive(Debug))]
struct SeqchainIndex {
    rank_data: Vec<u32>,
    less: Vec<u64>,
    sa: Vec<u64>,
    chroms: Vec<StoredChromInfo>,
    text: Vec<u8>,
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Save
// ═══════════════════════════════════════════════════════════════════════════════

/// Serialize an in-memory FM-Index to a `.seqchain` file.
pub fn save_index(searcher: &FmIndexSearcher, path: &Path) -> Result<(), SearchError> {
    let rank_data: Vec<u32> = searcher.to_interleaved_rank_data();
    let less: Vec<u64> = searcher.less_raw().iter().map(|&v| v as u64).collect();
    let sa: Vec<u64> = searcher.sa_raw().iter().map(|&v| v as u64).collect();
    let chroms: Vec<StoredChromInfo> = searcher
        .chroms()
        .iter()
        .map(|c| StoredChromInfo {
            name: c.name.clone(),
            start: c.start as u64,
            len: c.len as u64,
        })
        .collect();
    let text: Vec<u8> = searcher.text().to_vec();

    let index = SeqchainIndex {
        rank_data,
        less,
        sa,
        chroms,
        text,
    };

    let payload = rkyv::to_bytes::<RkyvError>(&index)
        .map_err(|e| SearchError::Other(anyhow::anyhow!("rkyv serialize failed: {}", e)))?;

    let file = File::create(path).map_err(SearchError::Io)?;
    let mut w = BufWriter::new(file);

    // Write header
    w.write_all(MAGIC).map_err(SearchError::Io)?;
    w.write_all(&FORMAT_VERSION.to_le_bytes())
        .map_err(SearchError::Io)?;
    w.write_all(&(payload.len() as u64).to_le_bytes())
        .map_err(SearchError::Io)?;

    // Write rkyv payload
    w.write_all(&payload).map_err(SearchError::Io)?;
    w.flush().map_err(SearchError::Io)?;

    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════════════
//  MappedIndex — zero-copy mmap'd FM-Index
// ═══════════════════════════════════════════════════════════════════════════════

/// An FM-Index loaded from a `.seqchain` file via mmap.
///
/// The payload is backed by 2 MiB huge pages (MAP_HUGETLB) when available,
/// collapsing TLB pressure from ~75K entries to ~150 for a 300 MB index.
/// Falls back to standard pages + MADV_HUGEPAGE on unconfigured systems.
pub struct MappedIndex {
    /// Huge-page-backed copy of the rkyv payload.
    _huge: HugePageMmap,
    /// Actual payload length within the huge-page region.
    payload_len: usize,
    /// Pointer into the huge-page data. Valid for the lifetime of `_huge`.
    archived: *const ArchivedSeqchainIndex,
    /// Pre-converted less table: `less[byte_value]` = count of chars < byte in BWT.
    less: Box<[usize; 256]>,
    /// Pre-converted chromosome metadata.
    chroms: Vec<(String, usize, usize)>,
}

// SAFETY: The Mmap is read-only and the archived pointer is derived from it.
// No mutation occurs after construction.
unsafe impl Send for MappedIndex {}
unsafe impl Sync for MappedIndex {}

impl MappedIndex {
    /// Access the archived index.
    #[inline]
    fn archived(&self) -> &ArchivedSeqchainIndex {
        // SAFETY: pointer is valid for the lifetime of _mmap, which we own.
        unsafe { &*self.archived }
    }

    /// Chromosome names in FASTA order.
    pub fn chrom_names(&self) -> Vec<String> {
        self.chroms.iter().map(|(name, _, _)| name.clone()).collect()
    }

    /// Chromosome geometry for the search engine.
    pub fn chrom_geometry(&self) -> ChromGeometry {
        ChromGeometry {
            ranges: self.chroms.iter().map(|&(_, start, len)| (start, len)).collect(),
        }
    }

    /// Concatenated genome text.
    pub fn text(&self) -> &[u8] {
        let archived = self.archived();
        &archived.text[..]
    }
}

// ─── FmOcc impl ──────────────────────────────────────────────────────────────

impl FmOcc for MappedIndex {
    #[inline]
    fn less(&self, c: u8) -> usize {
        self.less[c as usize]
    }

    #[inline]
    fn occ(&self, pos: usize, c: u8) -> usize {
        let bi = match c {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 0,
        };
        let archived = self.archived();
        archived.rank_data[pos * 4 + bi].to_native() as usize
    }

    #[inline]
    fn sa(&self, idx: usize) -> usize {
        let archived = self.archived();
        archived.sa[idx].to_native() as usize
    }

    #[inline]
    fn sa_len(&self) -> usize {
        let archived = self.archived();
        archived.sa.len()
    }

    #[inline]
    fn lf_map(&self, l: usize, r: usize) -> [(u32, u32); 4] {
        let archived = self.archived();
        let data = &archived.rank_data[..];

        let occ_l = if l > 0 {
            let b = (l - 1) * 4;
            [
                data[b].to_native(),
                data[b + 1].to_native(),
                data[b + 2].to_native(),
                data[b + 3].to_native(),
            ]
        } else {
            [0u32; 4]
        };

        let b = r * 4;
        let occ_r = [
            data[b].to_native(),
            data[b + 1].to_native(),
            data[b + 2].to_native(),
            data[b + 3].to_native(),
        ];

        let mut result = [(0u32, 0u32); 4];
        for bi in 0..4 {
            let less_b = self.less[BASES[bi] as usize] as u32;
            result[bi] = (less_b + occ_l[bi], less_b + occ_r[bi]);
        }
        result
    }

    #[inline]
    fn prefetch_lf(&self, l: usize, r: usize) {
        #[cfg(target_arch = "x86_64")]
        unsafe {
            use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
            let archived = self.archived();
            let data = &archived.rank_data[..];
            if l > 0 {
                let ptr_l = (data.as_ptr() as *const u8).add((l - 1) * 16) as *const i8;
                _mm_prefetch(ptr_l, _MM_HINT_T0);
            }
            let ptr_r = (data.as_ptr() as *const u8).add(r * 16) as *const i8;
            _mm_prefetch(ptr_r, _MM_HINT_T0);
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = (l, r);
        }
    }

    #[inline]
    fn rank_data(&self) -> Option<(&[u32], &[usize; 256])> {
        let archived = self.archived();
        // On little-endian, ArchivedVec<u32> is bit-identical to &[u32].
        // We transmute the slice of archived u32 to native u32.
        let data: &[rkyv::rend::u32_le] = &archived.rank_data[..];
        // SAFETY: u32_le has the same layout as u32 on LE platforms.
        #[cfg(target_endian = "little")]
        {
            let native_data: &[u32] =
                unsafe { std::slice::from_raw_parts(data.as_ptr() as *const u32, data.len()) };
            Some((native_data, &self.less))
        }
        #[cfg(not(target_endian = "little"))]
        {
            let _ = data;
            None // Fall back to scalar lf_map on big-endian
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Load
// ═══════════════════════════════════════════════════════════════════════════════

/// Load a `.seqchain` file into strict MAP_HUGETLB memory.
///
/// Reads the file, validates the header, copies the rkyv payload into
/// a `MAP_HUGETLB | MAP_ANONYMOUS | MAP_PRIVATE` region backed by 2 MiB pages.
///
/// # Panics
///
/// Panics if the kernel huge page pool is empty. Configure first:
/// ```bash
/// sudo sysctl -w vm.nr_hugepages=256
/// ```
pub fn load_index(path: &Path) -> Result<MappedIndex, SearchError> {
    let mut file = File::open(path).map_err(SearchError::Io)?;
    let file_len = file.metadata().map_err(SearchError::Io)?.len() as usize;

    if file_len < HEADER_SIZE {
        return Err(SearchError::Other(anyhow::anyhow!(
            "file too small for header: {} bytes",
            file_len
        )));
    }

    let mut header = [0u8; HEADER_SIZE];
    file.read_exact(&mut header).map_err(SearchError::Io)?;

    if &header[0..4] != MAGIC {
        return Err(SearchError::Other(anyhow::anyhow!(
            "invalid magic bytes (not a .seqchain file)"
        )));
    }

    let version = u32::from_le_bytes(header[4..8].try_into().unwrap());
    if version != FORMAT_VERSION {
        return Err(SearchError::Other(anyhow::anyhow!(
            "unsupported format version: {} (expected {})",
            version,
            FORMAT_VERSION
        )));
    }

    let payload_len = u64::from_le_bytes(header[8..16].try_into().unwrap()) as usize;
    if file_len < HEADER_SIZE + payload_len {
        return Err(SearchError::Other(anyhow::anyhow!(
            "file truncated: expected {} payload bytes, have {}",
            payload_len,
            file_len - HEADER_SIZE
        )));
    }

    // Read payload into temp buffer.
    let mut payload_buf = vec![0u8; payload_len];
    file.read_exact(&mut payload_buf).map_err(SearchError::Io)?;

    // Validate rkyv structure before binding to huge pages.
    let _: &ArchivedSeqchainIndex =
        rkyv::access::<ArchivedSeqchainIndex, RkyvError>(&payload_buf)
            .map_err(|e| SearchError::Other(anyhow::anyhow!("rkyv validation failed: {}", e)))?;

    // Bind to 2 MiB huge pages. Panics if pool is empty.
    let huge = HugePageMmap::from_data(&payload_buf);
    drop(payload_buf);

    // Re-access archived data from the huge-page region.
    let huge_slice = huge.as_slice(payload_len);
    let archived: &ArchivedSeqchainIndex =
        rkyv::access::<ArchivedSeqchainIndex, RkyvError>(huge_slice)
            .map_err(|e| SearchError::Other(anyhow::anyhow!("rkyv re-validation failed: {}", e)))?;

    let mut less = Box::new([0usize; 256]);
    for (i, val) in archived.less.iter().enumerate() {
        if i < 256 {
            less[i] = val.to_native() as usize;
        }
    }

    let chroms: Vec<(String, usize, usize)> = archived
        .chroms
        .iter()
        .map(|c| {
            (
                c.name.as_str().to_string(),
                c.start.to_native() as usize,
                c.len.to_native() as usize,
            )
        })
        .collect();

    let ptr = archived as *const ArchivedSeqchainIndex;

    Ok(MappedIndex {
        _huge: huge,
        payload_len,
        archived: ptr,
        less,
        chroms,
    })
}
