#!/usr/bin/env python3
"""Parquet sink benchmark — search_to_parquet vs search_batch.

Isolates the search + materialization path from the SeqChain pipeline.
Generates all guides up front, then benchmarks:
  1. search_batch  → Python lists (baseline)
  2. search_to_parquet → Arrow/Parquet (zero-allocation)
"""

import time
from pathlib import Path

from seqchain.io.genome import load_genbank
from seqchain.recipes import load_preset
from seqchain.recipes.crispr import interpret_guides
from seqchain.transform.regex import regex_map

from needletail import FmIndex

# ─── Paths ────────────────────────────────────────────────────────────

SEQCHAIN_ROOT = Path.home() / "Git" / "SeqChain"
GENOME_PATH = SEQCHAIN_ROOT / "tests" / "data" / "saccer3" / "sacCer3.gbff"
INDEX_DIR = Path(__file__).parent / "test" / "fixtures"
INDEX_PATH = str(INDEX_DIR / "sacCer3.seqchain")

MM = 2


def main():
    # ── Load genome + extract all guides ──────────────────────────────
    print("Loading genome + extracting guides …", flush=True)
    t0 = time.perf_counter()
    genome = load_genbank(GENOME_PATH)
    preset = load_preset("spcas9")

    hits = regex_map(genome.sequences, preset.pam, topologies=genome.topologies)
    guides = interpret_guides(hits, genome.sequences, preset, topologies=genome.topologies)

    # Materialize all guides into a list (so we can reuse them).
    all_guides = list(guides)
    queries = [g.tags["spacer"] for g in all_guides]
    t_prep = time.perf_counter() - t0
    print(f"  {len(queries):,} guides extracted in {t_prep:.2f}s")

    # ── Load index + warm seeds ───────────────────────────────────────
    print("\nLoading FM-Index …", flush=True)
    idx = FmIndex.load(INDEX_PATH)
    t0 = time.perf_counter()
    idx.warm_seeds(INDEX_PATH)
    t_warm = time.perf_counter() - t0
    print(f"  warm_seeds: {t_warm:.2f}s")

    # ── Benchmark 1: search_batch (Python lists) ─────────────────────
    print(f"\n{'═' * 60}")
    print(f"  BENCHMARK: search_batch → Python lists  (mm={MM})")
    print(f"{'═' * 60}")
    t0 = time.perf_counter()
    qi, pos, strand, scores = idx.search_batch(queries, mismatches=MM)
    t_batch = time.perf_counter() - t0
    print(f"  Hits: {len(qi):,}")
    print(f"  Time: {t_batch:.2f}s")

    # ── Benchmark 2: search_to_parquet (Arrow sink) ───────────────────
    print(f"\n{'═' * 60}")
    print(f"  BENCHMARK: search_to_parquet → Parquet  (mm={MM})")
    print(f"{'═' * 60}")
    out_path = "/tmp/needletail_bench.parquet"
    t0 = time.perf_counter()
    n_hits = idx.search_to_parquet(queries, out_path, mismatches=MM)
    t_parquet = time.perf_counter() - t0
    print(f"  Hits: {n_hits:,}")
    print(f"  Time: {t_parquet:.2f}s")

    parquet_size = Path(out_path).stat().st_size / (1024 * 1024)
    print(f"  File: {out_path} ({parquet_size:.1f} MiB)")

    # ── Delta ─────────────────────────────────────────────────────────
    print(f"\n{'─' * 60}")
    delta = (t_batch - t_parquet) / t_batch * 100
    print(f"  search_batch:      {t_batch:.2f}s")
    print(f"  search_to_parquet: {t_parquet:.2f}s")
    print(f"  Delta:             {delta:+.1f}%")
    print(f"{'─' * 60}")


if __name__ == "__main__":
    main()
