//! One-off pipeline: ZM56TN_1_WLP001 GBK + GFF → annotated Parquet.
//!
//! The GBK (AUGUSTUS predicted) has completely anonymous gene features.
//! The companion GFF has `transcript` lines with `GN=` gene symbols.
//! This script patches the genome's gene names from the GFF before running
//! the standard design_crispr_library pipeline.
//!
//! Matching key: (contig, start, end) in 0-based half-open coordinates.
//! GFF is 1-based closed → subtract 1 from start, keep end.
//!
//! Not a general GFF parser. Not a reusable tool. This specific file, once.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::Arc;

use needletail_core::engine::fm_index::ChromInfo;
use needletail_core::io::genbank::load_genbank;
use needletail_core::io::parquet_regions::ParquetFileSink;
use needletail_core::models::preset::{CRISPRPreset, FeatureConfig, NamingConfig};
use needletail_core::pipeline::design_crispr_library::{design_crispr_library, NullProgress};
use needletail_core::{build_seed_tiers, FmIndexSearcher, IndexHandle};

const GBK: &str = "/home/ryan.ward/Downloads/ZM56TN_1_WLP001_cells_1_reference.gbk";
const GFF: &str = "/home/ryan.ward/Downloads/ZM56TN_1_WLP001_cells_1_reference.gff";
const OUT: &str = "/home/ryan.ward/Downloads/ZM56TN_1_WLP001_cells_1_annotated.parquet";

fn main() {
    // ── 1. Parse GFF transcript lines → (contig, start0, end) → gene symbol ──
    //
    // GFF columns: contig · source · type · start(1-based) · end · score · strand · phase · attributes
    // We only care about `transcript` lines; gene lines are anonymous.
    // Attribute priority: GN= (gene symbol) → Name= (UniProt entry name)

    let mut gff_names: HashMap<(String, i64, i64), String> = HashMap::new();

    let gff_file = BufReader::new(File::open(GFF).expect("GFF not found"));
    for line in gff_file.lines() {
        let line = line.unwrap();
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.splitn(9, '\t').collect();
        if fields.len() < 9 || fields[2] != "transcript" {
            continue;
        }

        let contig = fields[0].to_string();
        let start0: i64 = fields[3].parse::<i64>().unwrap() - 1; // 1-based → 0-based
        let end: i64    = fields[4].parse::<i64>().unwrap();      // already exclusive in our model
        let attrs = fields[8];

        // Extract GN= first, fall back to Name=
        let name = extract_attr(attrs, "GN=")
            .or_else(|| extract_attr(attrs, "Name="))
            .unwrap_or_default();

        if !name.is_empty() {
            gff_names.insert((contig, start0, end), name);
        }
    }

    eprintln!("[gff] loaded {} named transcript entries", gff_names.len());

    // ── 2. Load GBK ──────────────────────────────────────────────────────────

    let naming = NamingConfig::default_config();
    let mut genome = load_genbank(Path::new(GBK), &naming.name_priority)
        .expect("Failed to load GBK");

    eprintln!("[gbk] {} contigs, {} features",
        genome.chromosomes.len(), genome.features.len());

    // ── 3. Patch gene names from GFF ─────────────────────────────────────────

    let mut patched = 0usize;
    for feat in &mut genome.features {
        if feat.name.is_empty() {
            let key = (feat.chrom.clone(), feat.start, feat.end);
            if let Some(name) = gff_names.get(&key) {
                feat.name = name.clone();
                patched += 1;
            }
        }
    }

    eprintln!("[patch] named {} / {} features", patched, genome.features.len());

    // ── 4. Build FM-Index ────────────────────────────────────────────────────

    let chroms: Vec<ChromInfo> = genome.chromosomes.iter()
        .map(|cm| ChromInfo { name: cm.name.clone(), start: cm.start, len: cm.len })
        .collect();
    let text = std::mem::take(&mut genome.text);
    let searcher = Arc::new(FmIndexSearcher::from_text(text, chroms).unwrap());

    let (ts, tl) = build_seed_tiers(&*searcher, searcher.text(), GBK).unwrap();
    let handle = IndexHandle::Built(searcher);

    // ── 5. Run pipeline ──────────────────────────────────────────────────────

    let preset = CRISPRPreset::by_name("spcas9").unwrap();
    let config = FeatureConfig::by_name("saccer3").unwrap();

    let mut sink = ParquetFileSink::create(Path::new(OUT)).unwrap();
    let result = design_crispr_library(
        &genome, &handle, Some(&ts), Some(&tl),
        &preset, &config, &NullProgress, &mut sink,
    ).unwrap();
    sink.finish().unwrap();

    eprintln!("[done] {} guides written → {}", result.guides_written, OUT);
}

/// Extract the value of a `KEY=value` pair from a GFF9 attributes string.
/// Stops at the next `;` or end of string.
fn extract_attr(attrs: &str, key: &str) -> Option<String> {
    let start = attrs.find(key)?;
    let value_start = start + key.len();
    let value_end = attrs[value_start..].find(';')
        .map(|i| value_start + i)
        .unwrap_or(attrs.len());
    let val = attrs[value_start..value_end].trim().to_string();
    if val.is_empty() { None } else { Some(val) }
}
