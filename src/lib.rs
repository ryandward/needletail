//! seqchain-native — PyO3 Python extension for width-first SIMD FM-Index search.
//!
//! Exposes `FmIndex` to Python. Construction is synchronous (builds the index
//! from a FASTA file). `search_batch` releases the GIL, runs the width-first
//! engine, and returns flat arrays directly to Python.

#![deny(clippy::all)]

mod error;
mod fm_index;
mod kmer_index;
mod simd_search;

use std::sync::Arc;

use bio::alphabets::dna;
use pyo3::prelude::*;

use crate::fm_index::FmIndexSearcher;
use crate::kmer_index::{KmerSeedTable, PosTable};
use crate::simd_search::{search_width_first, search_width_first_seeded, HitAccumulator};

// ─── FmIndex ─────────────────────────────────────────────────────────────────

/// FM-Index over a multi-FASTA genome.
///
/// Usage from Python:
/// ```python
/// from seqchain_native import FmIndex
///
/// idx = FmIndex("/path/to/genome.fa")
/// qi, pos, scores = idx.search_batch(["ATGATGATGATGATGATGATG"], mismatches=2)
/// ```
#[pyclass]
struct FmIndex {
    inner: Arc<FmIndexSearcher>,
    fasta_path: String,
    seed_table: Option<Arc<KmerSeedTable>>,
    pos_table: Option<Arc<PosTable>>,
}

#[pymethods]
impl FmIndex {
    /// Build an FM-Index from a FASTA file.
    ///
    /// Automatically builds (or loads from cache) the 10-mer seed table
    /// alongside the genome for accelerated seeding.
    #[new]
    fn new(fasta_path: &str) -> PyResult<Self> {
        let searcher = FmIndexSearcher::from_fasta(fasta_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?;
        let searcher = Arc::new(searcher);

        // Build or load the 10-mer seed table.
        let seed_table = KmerSeedTable::open_or_build(&*searcher, fasta_path, kmer_index::SEED_K)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                    "Failed to build/load k-mer seed table: {}",
                    e
                ))
            })?;

        // Build position table for verify-from-leading-seeds (pigeonhole).
        let pos_table = PosTable::build(searcher.text(), kmer_index::SEED_K);

        Ok(FmIndex {
            inner: searcher,
            fasta_path: fasta_path.to_string(),
            seed_table: Some(Arc::new(seed_table)),
            pos_table: Some(Arc::new(pos_table)),
        })
    }

    /// Search all queries against the genome with up to `mismatches` substitutions.
    ///
    /// All queries must be the same length.
    /// Releases the GIL for the duration of the search.
    ///
    /// Returns `(query_indices, positions, scores)` — three parallel lists.
    #[pyo3(signature = (queries, mismatches=0))]
    fn search_batch(
        &self,
        py: Python<'_>,
        queries: Vec<String>,
        mismatches: u8,
    ) -> PyResult<(Vec<u32>, Vec<u32>, Vec<f32>)> {
        if mismatches > 3 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "mismatches must be 0, 1, 2, or 3",
            ));
        }
        if queries.is_empty() {
            return Ok((vec![], vec![], vec![]));
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

        let chroms = self.inner.chrom_geometry();
        let index = self.inner.clone();
        let has_seeds = self.seed_table.is_some();

        // Use seeded path when seed table is available and queries are long enough.
        let has_pos = self.pos_table.is_some();
        if has_seeds && has_pos && query_len > kmer_index::SEED_K {
            let seed_table = self.seed_table.as_ref().unwrap().clone();
            let pos_table = self.pos_table.as_ref().unwrap().clone();

            let hits = py.allow_threads(move || {
                let genome_text = index.text();
                let mut acc = HitAccumulator::new();
                search_width_first_seeded(
                    &*index,
                    &*seed_table,
                    &*pos_table,
                    genome_text,
                    &queries_fwd,
                    &queries_rc,
                    query_len,
                    mismatches,
                    &chroms,
                    &mut acc,
                )
                .map(|()| acc)
            });

            match hits {
                Ok(acc) => Ok((acc.query_id, acc.position, acc.score)),
                Err(e) => Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                    e.to_string(),
                )),
            }
        } else {
            let hits = py.allow_threads(move || {
                let mut acc = HitAccumulator::new();
                search_width_first(
                    &*index,
                    &queries_fwd,
                    &queries_rc,
                    query_len,
                    mismatches,
                    &chroms,
                    &mut acc,
                );
                acc
            });

            Ok((hits.query_id, hits.position, hits.score))
        }
    }

    /// Exhaustive search (no seeding, no clog pruning) — for ground truth validation.
    #[pyo3(signature = (queries, mismatches=0))]
    fn search_batch_exhaustive(
        &self,
        py: Python<'_>,
        queries: Vec<String>,
        mismatches: u8,
    ) -> PyResult<(Vec<u32>, Vec<u32>, Vec<f32>)> {
        if mismatches > 3 || queries.is_empty() {
            return Ok((vec![], vec![], vec![]));
        }
        let query_len = queries[0].len();
        let queries_fwd: Vec<Vec<u8>> = queries
            .iter()
            .map(|q| q.bytes().map(|b| b.to_ascii_uppercase()).collect())
            .collect();
        let queries_rc: Vec<Vec<u8>> = queries_fwd.iter().map(|fwd| dna::revcomp(fwd)).collect();
        let chroms = self.inner.chrom_geometry();
        let index = self.inner.clone();
        let hits = py.allow_threads(move || {
            let mut acc = HitAccumulator::new();
            search_width_first(
                &*index, &queries_fwd, &queries_rc,
                query_len, mismatches, &chroms, &mut acc,
            );
            acc
        });
        Ok((hits.query_id, hits.position, hits.score))
    }

    /// Chromosome names in FASTA order (index = chrom_id).
    fn chrom_names(&self) -> Vec<String> {
        self.inner
            .chrom_names()
            .into_iter()
            .map(|s| s.to_string())
            .collect()
    }

    /// Chromosome geometry: list of `(start, length)` tuples.
    fn chrom_ranges(&self) -> Vec<(usize, usize)> {
        self.inner.chrom_geometry().ranges
    }
}

// ─── Module ──────────────────────────────────────────────────────────────────

#[pymodule]
fn seqchain_native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<FmIndex>()?;
    Ok(())
}
