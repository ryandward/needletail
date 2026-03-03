//! PyO3 bindings for Needletail — thin wrapper over needletail-core.
//!
//! This is the only crate that imports `pyo3`. It wraps the pure-Rust
//! `needletail_core` API into Python-accessible `#[pyclass]` types.

#![deny(clippy::all)]

use std::path::Path;
use std::sync::Arc;

use pyo3::prelude::*;

use needletail_core::{
    build_seed_tier_for_handle, build_seed_tiers, filter_hits_by_pam, prepare_queries,
    run_search_seeded, run_search_unseeded, select_tier, CompiledPam, FmIndexSearcher,
    GuideHits, IndexHandle, MappedIndex, PamDirection, SeedTier, SEED_K_LARGE, SEED_K_SMALL,
};
use needletail_core::operations::pam_scanner::{enrich_hits, find_pam_sites};

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
    inner: IndexHandle,
    fasta_path: Option<String>,
    tier_small: Option<SeedTier>,
    tier_large: Option<SeedTier>,
}

#[pymethods]
impl FmIndex {
    /// Build an FM-Index from a FASTA file.
    #[new]
    fn new(fasta_path: &str) -> PyResult<Self> {
        let searcher = FmIndexSearcher::from_fasta(fasta_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;
        let searcher = Arc::new(searcher);

        let (tier_small, tier_large) = build_seed_tiers(&*searcher, searcher.text(), fasta_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))?;

        Ok(FmIndex {
            inner: IndexHandle::Built(searcher),
            fasta_path: Some(fasta_path.to_string()),
            tier_small: Some(tier_small),
            tier_large: Some(tier_large),
        })
    }

    /// Build an FM-Index from a FASTA file and save it to a `.seqchain` file.
    #[staticmethod]
    fn build(fasta_path: &str, index_path: &str) -> PyResult<Self> {
        let searcher = FmIndexSearcher::from_fasta(fasta_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

        needletail_core::io::persist::save_index(&searcher, Path::new(index_path))
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;

        let searcher = Arc::new(searcher);

        let (tier_small, tier_large) = build_seed_tiers(&*searcher, searcher.text(), fasta_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))?;

        Ok(FmIndex {
            inner: IndexHandle::Built(searcher),
            fasta_path: Some(fasta_path.to_string()),
            tier_small: Some(tier_small),
            tier_large: Some(tier_large),
        })
    }

    /// Load a pre-built FM-Index from a `.seqchain` file via mmap.
    #[staticmethod]
    fn load(index_path: &str) -> PyResult<Self> {
        let mapped = needletail_core::io::persist::load_index(Path::new(index_path))
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;
        let mapped = Arc::new(mapped);

        Ok(FmIndex {
            inner: IndexHandle::Loaded(mapped),
            fasta_path: None,
            tier_small: None,
            tier_large: None,
        })
    }

    /// Build both seed tiers (K=10 and K=14) for seeded search acceleration.
    fn warm_seeds(&mut self, index_path: &str) -> PyResult<()> {
        let text = self.inner.text();

        if let Some(tier) = build_seed_tier_for_handle(&self.inner, text, index_path, SEED_K_SMALL)
        {
            self.tier_small = Some(tier);
        }

        if let Some(tier) = build_seed_tier_for_handle(&self.inner, text, index_path, SEED_K_LARGE)
        {
            self.tier_large = Some(tier);
        }

        Ok(())
    }

    /// Search all queries against the genome with up to `mismatches` substitutions.
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

        let (query_len, queries_fwd, queries_rc) = prepare_queries(&queries, mismatches)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))?;

        let chroms = self.inner.chrom_geometry();

        let compiled_pam = match pam {
            Some(p) => Some(
                CompiledPam::compile(p)
                    .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))?,
            ),
            None => None,
        };
        let direction = parse_direction(pam_direction)?;

        let mut acc = if let Some((tier, _k)) =
            select_tier(self.tier_small.as_ref(), self.tier_large.as_ref(), mismatches, query_len)
        {
            let seed_table = tier.seed_table.clone();
            let pos_table = tier.pos_table.clone();

            match &self.inner {
                IndexHandle::Built(index) => {
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
                IndexHandle::Loaded(index) => {
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
                IndexHandle::Built(index) => {
                    let index = index.clone();
                    py.allow_threads(move || {
                        run_search_unseeded(
                            &*index, &queries_fwd, &queries_rc,
                            query_len, mismatches, max_width, &chroms,
                        )
                    })
                }
                IndexHandle::Loaded(index) => {
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

        if let Some(ref compiled) = compiled_pam {
            let text = self.inner.text();
            let filter_chroms = self.inner.chrom_geometry();
            filter_hits_by_pam(
                &mut acc, text, &filter_chroms, compiled, direction,
                query_len, topologies.as_deref(),
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
        let (query_len, queries_fwd, queries_rc) = prepare_queries(&queries, mismatches)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))?;
        let chroms = self.inner.chrom_geometry();

        match &self.inner {
            IndexHandle::Built(index) => {
                let index = index.clone();
                let hits = py.allow_threads(move || {
                    run_search_unseeded(
                        &*index, &queries_fwd, &queries_rc,
                        query_len, mismatches, u32::MAX, &chroms,
                    )
                });
                Ok((hits.query_id, hits.position, hits.strand, hits.score))
            }
            IndexHandle::Loaded(index) => {
                let index = index.clone();
                let hits = py.allow_threads(move || {
                    run_search_unseeded(
                        &*index, &queries_fwd, &queries_rc,
                        query_len, mismatches, u32::MAX, &chroms,
                    )
                });
                Ok((hits.query_id, hits.position, hits.strand, hits.score))
            }
        }
    }

    /// Search and write results directly to a Parquet file.
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

        let (query_len, queries_fwd, queries_rc) = prepare_queries(&queries, mismatches)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))?;

        let chroms = self.inner.chrom_geometry();
        let out_path = std::path::PathBuf::from(output_path);

        let acc = if let Some((tier, _k)) =
            select_tier(self.tier_small.as_ref(), self.tier_large.as_ref(), mismatches, query_len)
        {
            let seed_table = tier.seed_table.clone();
            let pos_table = tier.pos_table.clone();

            match &self.inner {
                IndexHandle::Built(index) => {
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
                IndexHandle::Loaded(index) => {
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
                IndexHandle::Built(index) => {
                    let index = index.clone();
                    py.allow_threads(move || {
                        run_search_unseeded(
                            &*index, &queries_fwd, &queries_rc,
                            query_len, mismatches, max_width, &chroms,
                        )
                    })
                }
                IndexHandle::Loaded(index) => {
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

        let chroms_geo = self.inner.chrom_geometry();
        let n_rows = needletail_core::io::parquet_hits::hits_to_parquet(&acc, &chroms_geo, &out_path)
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
    #[pyo3(signature = (pam, spacer_len, pam_direction="downstream", topologies=None))]
    fn find_guides(
        &self,
        py: Python<'_>,
        pam: &str,
        spacer_len: usize,
        pam_direction: &str,
        topologies: Option<Vec<bool>>,
    ) -> PyResult<(Vec<u32>, Vec<u32>, Vec<i8>, Vec<u8>, Vec<u8>)> {
        let direction = parse_direction(pam_direction)?;

        if spacer_len == 0 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "spacer_len must be > 0",
            ));
        }

        let pam_owned = pam.to_string();
        let chroms = self.inner.chrom_geometry();

        let result = match &self.inner {
            IndexHandle::Built(index) => {
                let index = index.clone();
                py.allow_threads(move || {
                    let topo_ref = topologies.as_deref();
                    find_pam_sites(index.text(), &chroms, &pam_owned, spacer_len, direction, topo_ref)
                })
            }
            IndexHandle::Loaded(index) => {
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
    #[pyo3(signature = (pam, spacer_len, pam_direction="downstream", topologies=None))]
    fn scan_guides(
        &self,
        py: Python<'_>,
        pam: &str,
        spacer_len: usize,
        pam_direction: &str,
        topologies: Option<Vec<bool>>,
    ) -> PyResult<GuideBuffer> {
        let direction = parse_direction(pam_direction)?;

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
            IndexHandle::Built(index) => {
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
            IndexHandle::Loaded(index) => {
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

// ─── Helpers ─────────────────────────────────────────────────────────────────

fn parse_direction(s: &str) -> PyResult<PamDirection> {
    match s {
        "downstream" => Ok(PamDirection::Downstream),
        "upstream" => Ok(PamDirection::Upstream),
        _ => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "pam_direction must be 'downstream' or 'upstream'",
        )),
    }
}

// ─── Module ──────────────────────────────────────────────────────────────────

#[pymodule]
fn needletail(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<FmIndex>()?;
    m.add_class::<GuideBuffer>()?;
    Ok(())
}
