#!/usr/bin/env python3
"""Fully native SeqChain pipeline — Rust PAM scan + Rust PAM-validated scorer.

Composes the SeqChain generator chain with Rust-native guide scanning:

    scan_guides (Rust) → score_off_targets_fast → anchor_features
    → annotate_tracks → filter promoters → append_tss_distance
    → stream_json_export

All heavy-lifting in Rust: PAM scanning, spacer extraction, off-target search
(width-capped BWT), and PAM validation (IUPAC bitmask engine on mmap'd text).
Python is a thin membrane for Region construction and pipeline wiring.

Usage:
    python bench_annotate_saccer3.py -o sacCer3_native_idiomatic.json
    python bench_annotate_saccer3.py --mm 2 --max-width 16 -o output.json
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import Counter
from pathlib import Path

# ─── Zone 1: The atom ────────────────────────────────────────────────
from seqchain.region import Region

# ─── Zone 2: I/O ─────────────────────────────────────────────────────
from seqchain.io.genome import load_genbank
from seqchain.io.tracks import stream_json_export

# ─── Zone 4: Recipes ─────────────────────────────────────────────────
from seqchain.recipes import load_feature_config, load_preset
from seqchain.recipes.anchor import anchor_features
from seqchain.recipes.annotate import annotate_tracks, append_tss_distance

# ─── Native engine ───────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent))
from needletail import FmIndex
from python.needletail_guides import scan_guides, score_off_targets_fast


# ═════════════════════════════════════════════════════════════════════
#  Paths
# ═════════════════════════════════════════════════════════════════════

SEQCHAIN_ROOT = Path.home() / "Git" / "SeqChain"
GENOME_PATH = SEQCHAIN_ROOT / "tests" / "data" / "saccer3" / "sacCer3.gbff"

INDEX_DIR = Path(__file__).parent / "test" / "fixtures"
INDEX_PATH = str(INDEX_DIR / "sacCer3.seqchain")


def main():
    parser = argparse.ArgumentParser(
        description="SacCer3 SpCas9 annotation — native SIMD pipeline",
    )
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output JSON path (via stream_json_export)")
    parser.add_argument("--mm", type=int, default=2,
                        help="Mismatches for off-target search (default: 2)")
    parser.add_argument("--max-width", type=int, default=8,
                        help="BWT width cap per guide (default: 8)")
    args = parser.parse_args()

    t_total = time.perf_counter()

    # ═════════════════════════════════════════════════════════════════
    #  ZONE 2: LOAD
    # ═════════════════════════════════════════════════════════════════

    print("Loading SacCer3 genome …", flush=True)
    t0 = time.perf_counter()
    genome = load_genbank(GENOME_PATH)
    t_load = time.perf_counter() - t0

    total_bp = sum(genome.chrom_lengths.values())
    print(
        f"  {len(genome.sequences)} chromosomes · "
        f"{total_bp:,} bp · "
        f"{len(genome.genes()):,} genes · "
        f"loaded in {t_load:.2f}s"
    )

    # ═════════════════════════════════════════════════════════════════
    #  ZONE 4: CHEMISTRY + BIOLOGY
    # ═════════════════════════════════════════════════════════════════

    preset = load_preset("spcas9")
    config = load_feature_config("saccer3_features")
    print(f"Chemistry: {preset.name}  PAM={preset.pam}  spacer={preset.spacer_len}bp")
    print(f"Biology:   {config.organism}  ({len(config.features)} feature types)")

    # ═════════════════════════════════════════════════════════════════
    #  NATIVE ENGINE: mmap load + warm seeds
    # ═════════════════════════════════════════════════════════════════

    print(f"\nLoading native FM-Index …", flush=True)
    t0 = time.perf_counter()
    idx = FmIndex.load(INDEX_PATH)
    t_idx = time.perf_counter() - t0
    print(f"  mmap load: {t_idx*1e6:.0f} µs")

    print("Warming seed tables …", flush=True)
    t0 = time.perf_counter()
    idx.warm_seeds(INDEX_PATH)
    t_warm = time.perf_counter() - t0
    print(f"  warm_seeds: {t_warm:.2f}s")

    # ═════════════════════════════════════════════════════════════════
    #  PIPELINE ASSEMBLY — nothing computed yet
    # ═════════════════════════════════════════════════════════════════
    #
    #  scan_guides             → guide Regions (Rust: PAM scan + enrich)
    #    │
    #    ▼  score_off_targets_fast → scored Regions  (single Rust search_batch)
    #    ▼  annotate_tracks     → feature_type stamped      (sweep-line)
    #    ▼  filter promoters    → promoter guides only      (generator)
    #    ▼  append_tss_distance → signed_distance tag       (generator)
    #    ▼  stream_json_export  → DRAIN to disk
    #

    # Stage 1+2: Scan PAM sites + extract spacers (Rust, GIL released)
    guides = scan_guides(idx, preset, genome)

    # Stage 3: SIMD off-target scoring — single search_batch, PAM in Rust
    scored = score_off_targets_fast(
        guides, idx, preset.pam,
        spacer_len=preset.spacer_len,
        pam_direction=preset.pam_direction,
        mismatches=args.mm,
        topologies=genome.topologies,
        max_width=args.max_width,
    )

    # Stage 4: Tile features
    feature_tiles = anchor_features(
        genome.genes(), config, genome.chrom_lengths,
    )

    # Stage 5: Join guides × features via sweep-line
    annotated = annotate_tracks(scored, feature_tiles)

    # Stage 6: Filter promoters
    promoter_guides = (
        g for g in annotated
        if g.tags.get("feature_type") == "promoter"
    )

    # Stage 7: Append TSS distance
    with_distance = append_tss_distance(promoter_guides)

    # ═════════════════════════════════════════════════════════════════
    #  ZONE 2: THE DRAIN — computation starts here
    # ═════════════════════════════════════════════════════════════════

    total_regions = 0
    gene_counts: Counter[str] = Counter()
    off_target_scores: list[float] = []

    def spy(stream):
        nonlocal total_regions
        for r in stream:
            total_regions += 1
            fn = r.tags.get("feature_name")
            if fn:
                gene_counts[fn] += 1
            off_target_scores.append(r.score)
            yield r

    print(f"\nDraining pipeline → {args.output} …", flush=True)
    t0 = time.perf_counter()

    n = stream_json_export(args.output, spy(with_distance))

    t_drain = time.perf_counter() - t0
    t_wall = time.perf_counter() - t_total

    # ═════════════════════════════════════════════════════════════════
    #  SUMMARY
    # ═════════════════════════════════════════════════════════════════

    unique_genes = len(gene_counts)
    avg_off = sum(off_target_scores) / len(off_target_scores) if off_target_scores else 0
    perfect = sum(1 for s in off_target_scores if s == 0.0)
    pct_perfect = 100 * perfect / len(off_target_scores) if off_target_scores else 0

    print(f"\n{'─' * 60}")
    print(f"  SacCer3 SpCas9 CRISPRi promoter library — NATIVE SIMD")
    print(f"  mm={args.mm}  max_width={args.max_width}")
    print(f"{'─' * 60}")
    print(f"  Promoter-targeting guides     : {n:>12,}")
    print(f"  Unique genes targeted         : {unique_genes:>12,}")
    print(f"  Average off-target score      : {avg_off:>12.2f}")
    print(f"  Guides with score = 0.0       : {perfect:>12,}  ({pct_perfect:.1f}%)")
    print(f"{'─' * 60}")
    print(f"  Genome load  : {t_load:.2f}s")
    print(f"  mmap load    : {t_idx*1e6:.0f} µs")
    print(f"  warm_seeds   : {t_warm:.2f}s")
    print(f"  Pipeline drain: {t_drain:.2f}s")
    print(f"  Wall-clock    : {t_wall:.2f}s")
    print(f"{'─' * 60}")


if __name__ == "__main__":
    main()
