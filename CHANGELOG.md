# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-03-01

### Added

- **Rust-native PAM validation in `search_batch`** — new optional parameters
  `pam`, `pam_direction`, and `topologies` enable in-Rust IUPAC bitmask PAM
  filtering. Eliminates Python-side PAM regex validation entirely.
  28x scoring speedup (111 s → 3.95 s for 893K spacers on SacCer3).
- **`score_off_targets_fast()`** — streaming drop-in replacement for
  `score_off_targets_native()` and `score_off_targets()`. Macro-chunked
  (100K regions via `islice`), deduplicates spacers per chunk, issues one
  `search_batch` call per chunk with in-Rust PAM validation. O(chunk)
  memory, Axis 8 (Generator Purity) compliant. Accepts `**kwargs` for
  `ScorerFn` signature parity with existing SeqChain call sites.
- **`scan_guides()`** — Rust-native PAM scanning + guide enrichment via
  `FmIndex.scan_guides()`. Replaces `regex_map() + interpret_guides()`.
  All coordinate math, sequence assembly, and SHA-256 guide ID hashing
  done in Rust; Python only constructs `Region` objects.
- **`GuideBuffer` PyO3 class** — SoA container returned by `scan_guides()`,
  exposes `__len__` + `__getitem__` for zero-copy lazy iteration from Python.
- **`geometry.rs`** — pure 1D coordinate algebra (`normalize`,
  `interval_envelope`, `is_low_side`, `complement_base`, `fetch_sequence`).
  All `#[inline(always)]`, zero dependencies, branchless modular arithmetic.
- **`chemistry.rs`** — molecular rules leaf module (`PamDirection`,
  `CompiledPam`, `BASE_MASK`, `generate_guide_id`). No genome coordinates,
  no search structures, no Python.
- **`python/needletail_guides.py`** — thin Python bridge providing
  `scan_guides()` and `score_off_targets_fast()` as SeqChain-compatible
  generators.
- `sha2` dependency for deterministic guide ID hashing.

### Changed

- **5-pillar codebase reorganization** aligned with the SeqChain import
  hierarchy:
  - `src/geometry.rs` — pure 1D math (leaf, no deps)
  - `src/chemistry.rs` — molecular rules (leaf, no deps)
  - `src/engine/` — `fm_index.rs`, `kmer_index.rs`, `simd_search.rs`
  - `src/io/` — `persist.rs`, `parquet_sink.rs`
  - `src/operations/` — `pam_scanner.rs`
  - `src/lib.rs` — orchestration & FFI only
- `search_batch` signature extended with `max_width`, `pam`,
  `pam_direction`, `topologies` parameters (all optional, backward
  compatible).
- `bench_annotate_saccer3.py` updated to use `score_off_targets_fast`
  and `--max-width` parameter.

### Removed

- `src/primitives.rs` — split into `geometry.rs` (coordinates) and
  `chemistry.rs` (molecular rules).
- Flat `src/` layout — replaced by pillar-based module hierarchy.

## [0.1.0] - 2026-02-28

### Added

- Initial FM-Index implementation with BlockRank, AVX2 vertical automaton.
- Width-capped BWT search with tunable mismatch-branch pruning.
- Software prefetch pipeline for depth-2 pipelining.
- K=10 and K=14 seed tables for two-tier acceleration.
- rkyv zero-copy persistence with `MAP_HUGETLB` support.
- PyO3 bindings: `FmIndex.build()`, `.load()`, `.search_batch()`.
- Parquet output via Arrow zero-copy pipeline.
- Verified bit-identical results against Bowtie 1 at all mismatch levels.
