//! Orchestration & FFI — the only layer that knows Python exists.
//!
//! Holds the PyO3 `#[pyclass]` definitions (`FmIndex`, `GuideBuffer`),
//! manages the GIL (`py.allow_threads`), and orchestrates calls into the
//! five foundation pillars:
//!
//! ```text
//!   geometry    ← pure 1D coordinate math (leaf, no deps)
//!   chemistry   ← molecular rules: PAM masks, guide IDs (leaf, no deps)
//!   engine/     ← BWT, SIMD search, k-mer seeds (depends on geometry)
//!   io/         ← mmap persistence, Parquet sink (depends on engine)
//!   operations/ ← PAM scanning, off-target validation (depends on all above)
//! ```
//!
//! This file is the only place where `pyo3` is imported.

#![deny(clippy::all)]

// ─── The five pillars ────────────────────────────────────────────────────────
mod chemistry;
mod engine;
mod error;
mod geometry;
mod io;
mod operations;

use std::path::Path;
use std::sync::Arc;

use bio::alphabets::dna;
use pyo3::prelude::*;

use crate::chemistry::{CompiledPam, PamDirection};
use crate::engine::fm_index::FmIndexSearcher;
use crate::engine::kmer_index::{KmerSeedTable, PosTable, SEED_K_LARGE, SEED_K_SMALL};
use crate::engine::simd_search::{
    search_width_first, search_width_first_seeded, ChromGeometry, FmOcc, HitAccumulator,
};
use crate::io::persist::MappedIndex;
use crate::operations::pam_scanner::{
    enrich_hits, filter_hits_by_pam, find_pam_sites, GuideHits,
};

// ─── IndexInner ──────────────────────────────────────────────────────────────

enum IndexInner {
    Built(Arc<FmIndexSearcher>),
    Loaded(Arc<MappedIndex>),
}

impl IndexInner {
    fn chrom_geometry(&self) -> ChromGeometry {
        match self {
            IndexInner::Built(s) => s.chrom_geometry(),
            IndexInner::Loaded(m) => m.chrom_geometry(),
        }
    }

    fn chrom_names(&self) -> Vec<String> {
        match self {
            IndexInner::Built(s) => s.chrom_names().into_iter().map(|s| s.to_string()).collect(),
            IndexInner::Loaded(m) => m.chrom_names(),
        }
    }

    fn text(&self) -> &[u8] {
        match self {
            IndexInner::Built(s) => s.text(),
            IndexInner::Loaded(m) => m.text(),
        }
    }
}

// ─── Seed tier: paired KmerSeedTable + PosTable for a given K ────────────────

struct SeedTier {
    seed_table: Arc<KmerSeedTable>,
    pos_table: Arc<PosTable>,
}

// ─── Generic search helper ───────────────────────────────────────────────────

fn run_search_unseeded<I: FmOcc + Send + Sync>(
    index: &I,
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    query_len: usize,
    mismatches: u8,
    max_width: u32,
    chroms: &ChromGeometry,
) -> HitAccumulator {
    let mut acc = HitAccumulator::new();
    search_width_first(
        index,
        queries_fwd,
        queries_rc,
        query_len,
        mismatches,
        max_width,
        chroms,
        &mut acc,
    );
    acc
}

fn run_search_seeded<I: FmOcc + Send + Sync>(
    index: &I,
    seed_table: &KmerSeedTable,
    pos_table: &PosTable,
    genome_text: &[u8],
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    query_len: usize,
    mismatches: u8,
    max_width: u32,
    chroms: &ChromGeometry,
) -> Result<HitAccumulator, crate::error::SearchError> {
    use rayon::prelude::*;

    let n = queries_fwd.len();
    if n == 0 {
        return Ok(HitAccumulator::new());
    }

    // Split queries across rayon threads for parallel BWT processing.
    // Each chunk runs the full seeded pipeline independently.
    let n_threads = rayon::current_num_threads();
    let chunk_size = ((n + n_threads - 1) / n_threads).max(256);

    let chunks: Vec<(usize, usize)> = (0..n)
        .step_by(chunk_size)
        .map(|s| (s, (s + chunk_size).min(n)))
        .collect();

    let results: Result<Vec<HitAccumulator>, crate::error::SearchError> = chunks
        .par_iter()
        .map(|&(start, end)| {
            let mut acc = HitAccumulator::new();
            search_width_first_seeded(
                index,
                seed_table,
                pos_table,
                genome_text,
                &queries_fwd[start..end],
                &queries_rc[start..end],
                query_len,
                mismatches,
                max_width,
                chroms,
                &mut acc,
            )?;
            // Remap chunk-local query IDs to global IDs.
            for qi in &mut acc.query_id {
                *qi += start as u32;
            }
            Ok(acc)
        })
        .collect();

    let accs = results?;

    // Merge all chunk results.
    let total: usize = accs.iter().map(|a| a.query_id.len()).sum();
    let mut merged = HitAccumulator::new();
    merged.query_id.reserve(total);
    merged.position.reserve(total);
    merged.strand.reserve(total);
    merged.score.reserve(total);
    for acc in accs {
        merged.query_id.extend(acc.query_id);
        merged.position.extend(acc.position);
        merged.strand.extend(acc.strand);
        merged.score.extend(acc.score);
    }

    Ok(merged)
}

// ─── FmIndex ─────────────────────────────────────────────────────────────────

/// FM-Index over a multi-FASTA genome.
///
/// Usage from Python:
/// ```python
/// from needletail import FmIndex
///
/// idx = FmIndex("/path/to/genome.fa")              # build in-memory
/// idx = FmIndex.build("/path/to/genome.fa", "g.sc") # build + save
/// idx = FmIndex.load("g.sc")                        # mmap load, <50ms
/// qi, pos, strand, scores = idx.search_batch(["ATGATG"], mismatches=2)
/// ```
#[pyclass]
struct FmIndex {
    inner: IndexInner,
    fasta_path: Option<String>,
    /// K=10 tier: 8 MiB seed table, non-overlapping for L=20, fast at mm ≤ 2.
    tier_small: Option<SeedTier>,
    /// K=14 tier: 2.1 GiB seed table, deeper seed for mm = 3.
    tier_large: Option<SeedTier>,
}

impl FmIndex {
    /// Validate query batch and prepare forward + revcomp sequences.
    fn prepare_queries(
        queries: &[String],
        mismatches: u8,
    ) -> PyResult<(usize, Vec<Vec<u8>>, Vec<Vec<u8>>)> {
        if mismatches > 3 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "mismatches must be 0, 1, 2, or 3",
            ));
        }
        if queries.is_empty() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "queries must not be empty",
            ));
        }

        let query_len = queries[0].len();
        if query_len == 0 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "query length must be > 0",
            ));
        }
        for (i, q) in queries.iter().enumerate() {
            if q.len() != query_len {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "query {} has length {} but expected {}",
                    i,
                    q.len(),
                    query_len,
                )));
            }
        }

        let queries_fwd: Vec<Vec<u8>> = queries
            .iter()
            .map(|q| q.bytes().map(|b| b.to_ascii_uppercase()).collect())
            .collect();
        let queries_rc: Vec<Vec<u8>> = queries_fwd.iter().map(|fwd| dna::revcomp(fwd)).collect();

        Ok((query_len, queries_fwd, queries_rc))
    }

    /// Build both seed tiers (K=10 and K=14) from an FM-Index.
    fn build_tiers<I: FmOcc>(index: &I, text: &[u8], fasta_path: &str) -> PyResult<(SeedTier, SeedTier)> {
        // K=10 tier (8 MiB)
        let seed_small = KmerSeedTable::open_or_build(index, fasta_path, SEED_K_SMALL)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                format!("Failed to build/load K={} seed table: {}", SEED_K_SMALL, e),
            ))?;
        let pos_small = PosTable::build(text, SEED_K_SMALL);

        // K=14 tier (2.1 GiB)
        let seed_large = KmerSeedTable::open_or_build(index, fasta_path, SEED_K_LARGE)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                format!("Failed to build/load K={} seed table: {}", SEED_K_LARGE, e),
            ))?;
        let pos_large = PosTable::build(text, SEED_K_LARGE);

        Ok((
            SeedTier { seed_table: Arc::new(seed_small), pos_table: Arc::new(pos_small) },
            SeedTier { seed_table: Arc::new(seed_large), pos_table: Arc::new(pos_large) },
        ))
    }

    /// Select the appropriate seed tier based on mismatch level and query length.
    ///
    /// Dynamic routing:
    ///   mm ≤ 2 → K=10 (non-overlapping, seed_mm=1, 43 variants, 10 BWT steps)
    ///   mm = 3 → K=14 ONLY if 2K ≤ L (non-overlapping), else K=10
    ///
    /// Pigeonhole coverage requires non-overlapping segments. For K=14, L=20:
    /// segments overlap by 8bp, so 3 mismatches can all land in the overlap,
    /// making both seeds have 3mm (exceeding seed_mm=2). K=10 is safe because
    /// 2×10 ≤ 20, giving non-overlapping coverage with seed_mm=1.
    fn select_tier(&self, mismatches: u8, query_len: usize) -> Option<(&SeedTier, usize)> {
        if mismatches <= 2 {
            // Prefer K=10: non-overlapping segments for L=20, fast seeding
            if let Some(ref tier) = self.tier_small {
                if query_len > SEED_K_SMALL {
                    return Some((tier, SEED_K_SMALL));
                }
            }
        } else {
            // mm=3: use K=14 ONLY if segments are non-overlapping (2K ≤ L).
            // Otherwise pigeonhole coverage fails — fall through to K=10.
            if let Some(ref tier) = self.tier_large {
                if query_len >= 2 * SEED_K_LARGE {
                    return Some((tier, SEED_K_LARGE));
                }
            }
        }
        // Fallback: try whichever tier fits the query
        if let Some(ref tier) = self.tier_small {
            if query_len > SEED_K_SMALL {
                return Some((tier, SEED_K_SMALL));
            }
        }
        None
    }
}

#[pymethods]
impl FmIndex {
    /// Build an FM-Index from a FASTA file.
    #[new]
    fn new(fasta_path: &str) -> PyResult<Self> {
        let searcher = FmIndexSearcher::from_fasta(fasta_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;
        let searcher = Arc::new(searcher);

        let (tier_small, tier_large) = Self::build_tiers(&*searcher, searcher.text(), fasta_path)?;

        Ok(FmIndex {
            inner: IndexInner::Built(searcher),
            fasta_path: Some(fasta_path.to_string()),
            tier_small: Some(tier_small),
            tier_large: Some(tier_large),
        })
    }

    /// Build an FM-Index from a FASTA file and save it to a `.seqchain` file.
    ///
    /// Returns the index ready for searching.
    #[staticmethod]
    fn build(fasta_path: &str, index_path: &str) -> PyResult<Self> {
        let searcher = FmIndexSearcher::from_fasta(fasta_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

        io::persist::save_index(&searcher, Path::new(index_path))
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

        let searcher = Arc::new(searcher);

        let (tier_small, tier_large) = Self::build_tiers(&*searcher, searcher.text(), fasta_path)?;

        Ok(FmIndex {
            inner: IndexInner::Built(searcher),
            fasta_path: Some(fasta_path.to_string()),
            tier_small: Some(tier_small),
            tier_large: Some(tier_large),
        })
    }

    /// Load a pre-built FM-Index from a `.seqchain` file via mmap.
    ///
    /// This is the fast path: <50ms cold start for SacCer3.
    /// KmerSeedTable and PosTable are NOT rebuilt on load — searches use the
    /// unseeded path, which is still fast (in-process, no subprocess overhead).
    /// Call `warm_seeds()` after load to enable seeded search if needed.
    #[staticmethod]
    fn load(index_path: &str) -> PyResult<Self> {
        let mapped = io::persist::load_index(Path::new(index_path))
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;
        let mapped = Arc::new(mapped);

        Ok(FmIndex {
            inner: IndexInner::Loaded(mapped),
            fasta_path: None,
            tier_small: None,
            tier_large: None,
        })
    }

    /// Build both seed tiers (K=10 and K=14) for seeded search acceleration.
    ///
    /// Call this after `load()` to enable the seeded search path.
    /// Builds K=10 (8 MiB, instant) and K=14 (2.1 GiB, ~3s) seed tables,
    /// plus K=10 and K=14 PosTables from the mmap'd genome text.
    fn warm_seeds(&mut self, index_path: &str) -> PyResult<()> {
        let text = self.inner.text();
        let idx_path = Path::new(index_path);
        let stem = idx_path
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy();
        let parent = idx_path.parent().unwrap_or_else(|| Path::new("."));

        // Helper: build or load a KmerSeedTable for a given K
        let build_seed_table = |k: usize| -> Option<Arc<KmerSeedTable>> {
            let kmer_path = parent.join(format!("{}.{}mer.idx", stem, k));
            let table = match &self.inner {
                IndexInner::Built(s) => {
                    KmerSeedTable::open_or_build(&**s, "<mmap>", k).ok()
                }
                IndexInner::Loaded(m) => {
                    if kmer_path.exists() {
                        KmerSeedTable::open(&kmer_path, k).ok()
                    } else {
                        engine::kmer_index::build_kmer_index(&**m, k, &kmer_path)
                            .ok()
                            .and_then(|()| KmerSeedTable::open(&kmer_path, k).ok())
                    }
                }
            };
            table.map(Arc::new)
        };

        // K=10 tier
        if let Some(seed_small) = build_seed_table(SEED_K_SMALL) {
            let pos_small = Arc::new(PosTable::build(text, SEED_K_SMALL));
            self.tier_small = Some(SeedTier { seed_table: seed_small, pos_table: pos_small });
        }

        // K=14 tier
        if let Some(seed_large) = build_seed_table(SEED_K_LARGE) {
            let pos_large = Arc::new(PosTable::build(text, SEED_K_LARGE));
            self.tier_large = Some(SeedTier { seed_table: seed_large, pos_table: pos_large });
        }

        Ok(())
    }

    /// Search all queries against the genome with up to `mismatches` substitutions.
    ///
    /// All queries must be the same length.
    /// Releases the GIL for the duration of the search.
    ///
    /// Dynamic routing:
    ///   mm ≤ 2 → K=10 seed tier (non-overlapping, fast seeding)
    ///   mm = 3 → K=14 seed tier (deeper seed, fewer BWT steps)
    ///
    /// When `pam` is provided, hits are filtered in Rust using the IUPAC bitmask
    /// engine: only positions with a valid adjacent PAM survive. This replaces
    /// Python-side PAM regex validation entirely.
    ///
    /// Returns `(query_indices, positions, strands, scores)` — four parallel lists.
    #[pyo3(signature = (queries, mismatches=0, max_width=u32::MAX, pam=None, pam_direction="downstream", topologies=None))]
    fn search_batch(
        &self,
        py: Python<'_>,
        queries: Vec<String>,
        mismatches: u8,
        max_width: u32,
        pam: Option<&str>,
        pam_direction: &str,
        topologies: Option<Vec<bool>>,
    ) -> PyResult<(Vec<u32>, Vec<u32>, Vec<bool>, Vec<f32>)> {
        if queries.is_empty() {
            return Ok((vec![], vec![], vec![], vec![]));
        }

        let (query_len, queries_fwd, queries_rc) =
            Self::prepare_queries(&queries, mismatches)?;

        let chroms = self.inner.chrom_geometry();

        // Compile PAM filter if provided.
        let compiled_pam = match pam {
            Some(p) => Some(
                CompiledPam::compile(p)
                    .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))?,
            ),
            None => None,
        };
        let direction = match pam_direction {
            "downstream" => PamDirection::Downstream,
            "upstream" => PamDirection::Upstream,
            _ => {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    "pam_direction must be 'downstream' or 'upstream'",
                ));
            }
        };

        // ── Search (GIL released) ────────────────────────────────────
        let mut acc = if let Some((tier, _k)) = self.select_tier(mismatches, query_len) {
            let seed_table = tier.seed_table.clone();
            let pos_table = tier.pos_table.clone();

            match &self.inner {
                IndexInner::Built(index) => {
                    let index = index.clone();
                    py.allow_threads(move || {
                        run_search_seeded(
                            &*index,
                            &seed_table,
                            &pos_table,
                            index.text(),
                            &queries_fwd,
                            &queries_rc,
                            query_len,
                            mismatches,
                            max_width,
                            &chroms,
                        )
                    })
                    .map_err(|e| {
                        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string())
                    })?
                }
                IndexInner::Loaded(index) => {
                    let index = index.clone();
                    py.allow_threads(move || {
                        run_search_seeded(
                            &*index,
                            &seed_table,
                            &pos_table,
                            index.text(),
                            &queries_fwd,
                            &queries_rc,
                            query_len,
                            mismatches,
                            max_width,
                            &chroms,
                        )
                    })
                    .map_err(|e| {
                        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string())
                    })?
                }
            }
        } else {
            // ── Unseeded fallback ────────────────────────────────────
            match &self.inner {
                IndexInner::Built(index) => {
                    let index = index.clone();
                    py.allow_threads(move || {
                        run_search_unseeded(
                            &*index,
                            &queries_fwd,
                            &queries_rc,
                            query_len,
                            mismatches,
                            max_width,
                            &chroms,
                        )
                    })
                }
                IndexInner::Loaded(index) => {
                    let index = index.clone();
                    py.allow_threads(move || {
                        run_search_unseeded(
                            &*index,
                            &queries_fwd,
                            &queries_rc,
                            query_len,
                            mismatches,
                            max_width,
                            &chroms,
                        )
                    })
                }
            }
        };

        // ── PAM filter (pure Rust, runs after search) ────────────────
        if let Some(ref compiled) = compiled_pam {
            let text = self.inner.text();
            let filter_chroms = self.inner.chrom_geometry();
            filter_hits_by_pam(
                &mut acc,
                text,
                &filter_chroms,
                compiled,
                direction,
                query_len,
                topologies.as_deref(),
            );
        }

        Ok((acc.query_id, acc.position, acc.strand, acc.score))
    }

    /// Exhaustive search (no seeding, no clog pruning) — for ground truth validation.
    #[pyo3(signature = (queries, mismatches=0))]
    fn search_batch_exhaustive(
        &self,
        py: Python<'_>,
        queries: Vec<String>,
        mismatches: u8,
    ) -> PyResult<(Vec<u32>, Vec<u32>, Vec<bool>, Vec<f32>)> {
        if mismatches > 3 || queries.is_empty() {
            return Ok((vec![], vec![], vec![], vec![]));
        }
        let query_len = queries[0].len();
        let queries_fwd: Vec<Vec<u8>> = queries
            .iter()
            .map(|q| q.bytes().map(|b| b.to_ascii_uppercase()).collect())
            .collect();
        let queries_rc: Vec<Vec<u8>> = queries_fwd.iter().map(|fwd| dna::revcomp(fwd)).collect();
        let chroms = self.inner.chrom_geometry();

        match &self.inner {
            IndexInner::Built(index) => {
                let index = index.clone();
                let hits = py.allow_threads(move || {
                    run_search_unseeded(
                        &*index,
                        &queries_fwd,
                        &queries_rc,
                        query_len,
                        mismatches,
                        u32::MAX,
                        &chroms,
                    )
                });
                Ok((hits.query_id, hits.position, hits.strand, hits.score))
            }
            IndexInner::Loaded(index) => {
                let index = index.clone();
                let hits = py.allow_threads(move || {
                    run_search_unseeded(
                        &*index,
                        &queries_fwd,
                        &queries_rc,
                        query_len,
                        mismatches,
                        u32::MAX,
                        &chroms,
                    )
                });
                Ok((hits.query_id, hits.position, hits.strand, hits.score))
            }
        }
    }

    /// Search and write results directly to a Parquet file.
    ///
    /// Zero-allocation pipeline: wavefront survivors push raw bits into Arrow
    /// columnar builders, then flush to a `.parquet` sidecar in one shot.
    /// Returns the number of hits written.
    #[pyo3(signature = (queries, output_path, mismatches=0, max_width=u32::MAX))]
    fn search_to_parquet(
        &self,
        py: Python<'_>,
        queries: Vec<String>,
        output_path: &str,
        mismatches: u8,
        max_width: u32,
    ) -> PyResult<usize> {
        if queries.is_empty() {
            return Ok(0);
        }

        let (query_len, queries_fwd, queries_rc) =
            Self::prepare_queries(&queries, mismatches)?;

        let chroms = self.inner.chrom_geometry();
        let out_path = std::path::PathBuf::from(output_path);

        // Run search (reuse existing path to get HitAccumulator).
        let acc = if let Some((tier, _k)) = self.select_tier(mismatches, query_len) {
            let seed_table = tier.seed_table.clone();
            let pos_table = tier.pos_table.clone();

            match &self.inner {
                IndexInner::Built(index) => {
                    let index = index.clone();
                    py.allow_threads(move || {
                        run_search_seeded(
                            &*index, &seed_table, &pos_table, index.text(),
                            &queries_fwd, &queries_rc, query_len, mismatches,
                            max_width, &chroms,
                        )
                    })
                    .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?
                }
                IndexInner::Loaded(index) => {
                    let index = index.clone();
                    py.allow_threads(move || {
                        run_search_seeded(
                            &*index, &seed_table, &pos_table, index.text(),
                            &queries_fwd, &queries_rc, query_len, mismatches,
                            max_width, &chroms,
                        )
                    })
                    .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?
                }
            }
        } else {
            match &self.inner {
                IndexInner::Built(index) => {
                    let index = index.clone();
                    py.allow_threads(move || {
                        run_search_unseeded(
                            &*index, &queries_fwd, &queries_rc,
                            query_len, mismatches, max_width, &chroms,
                        )
                    })
                }
                IndexInner::Loaded(index) => {
                    let index = index.clone();
                    py.allow_threads(move || {
                        run_search_unseeded(
                            &*index, &queries_fwd, &queries_rc,
                            query_len, mismatches, max_width, &chroms,
                        )
                    })
                }
            }
        };

        // Sink directly to Parquet via Arrow builders — zero string allocation.
        let chroms_geo = self.inner.chrom_geometry();
        let n_rows = io::parquet_sink::hits_to_parquet(&acc, &chroms_geo, &out_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

        Ok(n_rows)
    }

    /// Chromosome names in FASTA order (index = chrom_id).
    fn chrom_names(&self) -> Vec<String> {
        self.inner.chrom_names()
    }

    /// Chromosome geometry: list of `(start, length)` tuples.
    fn chrom_ranges(&self) -> Vec<(usize, usize)> {
        self.inner.chrom_geometry().ranges
    }

    /// Scan the genome for CRISPR guide sites defined by a PAM motif.
    ///
    /// Returns `(chrom_ids, positions, strands, spacers, pam_seqs)` — five
    /// flat arrays. `spacers` and `pam_seqs` are packed byte buffers with
    /// stride `spacer_len` and `len(pam)` respectively.
    ///
    /// ```python
    /// c, p, s, spacers, pams = idx.find_guides("NGG", 20, "downstream")
    /// # spacers[i*20:(i+1)*20] is the i-th spacer sequence
    /// ```
    #[pyo3(signature = (pam, spacer_len, pam_direction="downstream", topologies=None))]
    fn find_guides(
        &self,
        py: Python<'_>,
        pam: &str,
        spacer_len: usize,
        pam_direction: &str,
        topologies: Option<Vec<bool>>,
    ) -> PyResult<(Vec<u32>, Vec<u32>, Vec<i8>, Vec<u8>, Vec<u8>)> {
        let direction = match pam_direction {
            "downstream" => PamDirection::Downstream,
            "upstream" => PamDirection::Upstream,
            _ => {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    "pam_direction must be 'downstream' or 'upstream'",
                ));
            }
        };

        if spacer_len == 0 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "spacer_len must be > 0",
            ));
        }

        let pam_owned = pam.to_string();
        let chroms = self.inner.chrom_geometry();

        let result = match &self.inner {
            IndexInner::Built(index) => {
                let index = index.clone();
                py.allow_threads(move || {
                    let topo_ref = topologies.as_deref();
                    find_pam_sites(index.text(), &chroms, &pam_owned, spacer_len, direction, topo_ref)
                })
            }
            IndexInner::Loaded(index) => {
                let index = index.clone();
                py.allow_threads(move || {
                    let topo_ref = topologies.as_deref();
                    find_pam_sites(index.text(), &chroms, &pam_owned, spacer_len, direction, topo_ref)
                })
            }
        };

        match result {
            Ok(hits) => Ok((hits.chrom_id, hits.position, hits.strand, hits.spacers, hits.pam_seqs)),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e)),
        }
    }

    /// Scan for CRISPR guide sites and return an enriched GuideBuffer.
    ///
    /// Performs PAM scanning and guide enrichment (footprints, guide_seqs,
    /// guide_ids) entirely in Rust with the GIL released. Returns a
    /// `GuideBuffer` that supports `len()` and `buf[i]` access from Python.
    #[pyo3(signature = (pam, spacer_len, pam_direction="downstream", topologies=None))]
    fn scan_guides(
        &self,
        py: Python<'_>,
        pam: &str,
        spacer_len: usize,
        pam_direction: &str,
        topologies: Option<Vec<bool>>,
    ) -> PyResult<GuideBuffer> {
        let direction = match pam_direction {
            "downstream" => PamDirection::Downstream,
            "upstream" => PamDirection::Upstream,
            _ => {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    "pam_direction must be 'downstream' or 'upstream'",
                ));
            }
        };

        if spacer_len == 0 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "spacer_len must be > 0",
            ));
        }

        let pam_owned = pam.to_string();
        let chroms = self.inner.chrom_geometry();
        let chrom_names = self.inner.chrom_names();
        let chrom_ranges = chroms.ranges.clone();

        let result: Result<GuideHits, String> = match &self.inner {
            IndexInner::Built(index) => {
                let index = index.clone();
                let names = chrom_names.clone();
                py.allow_threads(move || {
                    let topo_ref = topologies.as_deref();
                    let hits = find_pam_sites(
                        index.text(), &chroms, &pam_owned, spacer_len, direction, topo_ref,
                    )?;
                    Ok(enrich_hits(hits, &names, direction, topo_ref, &chrom_ranges))
                })
            }
            IndexInner::Loaded(index) => {
                let index = index.clone();
                let names = chrom_names.clone();
                py.allow_threads(move || {
                    let topo_ref = topologies.as_deref();
                    let hits = find_pam_sites(
                        index.text(), &chroms, &pam_owned, spacer_len, direction, topo_ref,
                    )?;
                    Ok(enrich_hits(hits, &names, direction, topo_ref, &chrom_ranges))
                })
            }
        };

        match result {
            Ok(guide_hits) => Ok(GuideBuffer { hits: guide_hits, chrom_names }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e)),
        }
    }
}

// ─── GuideBuffer ─────────────────────────────────────────────────────────────

#[pyclass]
struct GuideBuffer {
    hits: GuideHits,
    chrom_names: Vec<String>,
}

#[pymethods]
impl GuideBuffer {
    #[getter]
    fn count(&self) -> usize {
        self.hits.count
    }

    #[getter]
    fn spacer_len(&self) -> usize {
        self.hits.spacer_len
    }

    #[getter]
    fn pam_len(&self) -> usize {
        self.hits.pam_len
    }

    #[getter]
    fn chrom_names(&self) -> Vec<String> {
        self.chrom_names.clone()
    }

    fn __len__(&self) -> usize {
        self.hits.count
    }

    fn __getitem__(&self, py: Python<'_>, idx: isize) -> PyResult<PyObject> {
        let n = self.hits.count as isize;
        let i = if idx < 0 { idx + n } else { idx };
        if i < 0 || i >= n {
            return Err(PyErr::new::<pyo3::exceptions::PyIndexError, _>(
                "index out of range",
            ));
        }
        let i = i as usize;
        let sl = self.hits.spacer_len;
        let pl = self.hits.pam_len;
        let gl = sl + pl;

        let cid = self.hits.chrom_ids[i];
        let pam_pos = self.hits.pam_positions[i];
        let gs = self.hits.guide_starts[i];
        let ge = self.hits.guide_ends[i];
        let strand_str = if self.hits.strands[i] > 0 { "+" } else { "-" };
        let spacer = &self.hits.spacers[i * sl..(i + 1) * sl];
        let pam_seq = &self.hits.pam_seqs[i * pl..(i + 1) * pl];
        let guide_seq = &self.hits.guide_seqs[i * gl..(i + 1) * gl];
        let guide_id = std::str::from_utf8(&self.hits.guide_ids[i * 8..(i + 1) * 8])
            .unwrap_or("");

        Ok((
            cid, pam_pos, gs, ge,
            strand_str, spacer, pam_seq, guide_seq, guide_id,
        )
            .into_pyobject(py)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?
            .into_any()
            .unbind())
    }
}

// ─── Module ──────────────────────────────────────────────────────────────────

#[pymodule]
fn needletail(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<FmIndex>()?;
    m.add_class::<GuideBuffer>()?;
    Ok(())
}
