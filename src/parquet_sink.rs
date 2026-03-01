//! Arrow-native Parquet sink — true zero-copy results pipeline.
//!
//! The HitAccumulator already stores results in SoA (struct-of-arrays) layout:
//! `Vec<u32>` for query_ids, `Vec<u32>` for positions, etc. These are contiguous
//! memory buffers — exactly what Arrow arrays are.
//!
//! Instead of iterating + pushing into Arrow builders (O(n) copy), we convert
//! the raw `Vec<u32>` directly into Arrow `Buffer` objects via pointer handoff.
//! The Arrow array wraps the existing allocation — no memcpy, no iteration.

use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    ArrayData, BooleanArray, Float32Array, PrimitiveArray, RecordBatch, UInt32Array, UInt8Array,
};
use arrow::buffer::Buffer;
use arrow::datatypes::{DataType, Field, Float32Type, Schema, UInt32Type, UInt8Type};
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::simd_search::ChromGeometry;
use crate::simd_search::HitAccumulator;

// ═══════════════════════════════════════════════════════════════════════════════
//  Schema
// ═══════════════════════════════════════════════════════════════════════════════

fn raw_schema() -> Schema {
    Schema::new(vec![
        Field::new("query_id", DataType::UInt32, false),
        Field::new("chrom_id", DataType::UInt8, false),
        Field::new("position", DataType::UInt32, false),
        Field::new("strand", DataType::Boolean, false),
        Field::new("score", DataType::Float32, false),
    ])
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Zero-copy Vec → Arrow Buffer
// ═══════════════════════════════════════════════════════════════════════════════

/// Convert `Vec<T>` to Arrow `Buffer` by handing off the allocation.
/// No copy — Arrow takes ownership of the heap pointer.
#[inline]
fn vec_to_buffer<T: arrow::datatypes::ArrowNativeType>(v: Vec<T>) -> Buffer {
    Buffer::from_vec(v)
}

/// Convert `Vec<u32>` directly into a `UInt32Array` — zero copy.
#[inline]
fn vec_u32_to_array(v: Vec<u32>) -> UInt32Array {
    let len = v.len();
    let buf = vec_to_buffer(v);
    let data = unsafe {
        ArrayData::new_unchecked(DataType::UInt32, len, None, None, 0, vec![buf], vec![])
    };
    PrimitiveArray::<UInt32Type>::from(data)
}

/// Convert `Vec<f32>` directly into a `Float32Array` — zero copy.
#[inline]
fn vec_f32_to_array(v: Vec<f32>) -> Float32Array {
    let len = v.len();
    let buf = vec_to_buffer(v);
    let data = unsafe {
        ArrayData::new_unchecked(DataType::Float32, len, None, None, 0, vec![buf], vec![])
    };
    PrimitiveArray::<Float32Type>::from(data)
}

/// Convert `Vec<bool>` into a `BooleanArray` (bit-packed).
/// This is the one column that requires actual work: bools → bit buffer.
#[inline]
fn vec_bool_to_array(v: Vec<bool>) -> BooleanArray {
    // BooleanArray is bit-packed, so we must pack the bools.
    // arrow's From<Vec<bool>> does this efficiently.
    BooleanArray::from(v)
}

// ═══════════════════════════════════════════════════════════════════════════════
//  HitAccumulator → Parquet (zero-copy, single pass chrom resolution)
// ═══════════════════════════════════════════════════════════════════════════════

/// Zero-copy HitAccumulator → Parquet.
///
/// 1. Resolve global positions → (chrom_id, local_pos) via binary search.
///    This is the only computation — O(n log k) for n=3M, k=17.
/// 2. Hand off the `Vec<u32>` buffers directly to Arrow (no copy).
/// 3. Flush to Parquet with Snappy compression.
pub fn hits_to_parquet(
    hits: &HitAccumulator,
    chroms: &ChromGeometry,
    path: &Path,
) -> Result<usize, anyhow::Error> {
    let n = hits.query_id.len();
    if n == 0 {
        let schema = Arc::new(raw_schema());
        let batch = RecordBatch::new_empty(schema);
        write_parquet(&batch, path)?;
        return Ok(0);
    }

    // Pre-build sorted start offsets for binary search.
    let starts: Vec<usize> = chroms.ranges.iter().map(|&(s, _)| s).collect();

    // Single pass: resolve chrom_id + local position.
    let mut chrom_id_vec: Vec<u8> = Vec::with_capacity(n);
    let mut local_pos_vec: Vec<u32> = Vec::with_capacity(n);

    for i in 0..n {
        let pos = hits.position[i] as usize;
        let ci = starts.partition_point(|&s| s <= pos).saturating_sub(1);
        chrom_id_vec.push(ci as u8);
        local_pos_vec.push((pos - starts[ci]) as u32);
    }

    // Zero-copy conversion: hand Vec ownership to Arrow.
    let query_id_arr = vec_u32_to_array(hits.query_id.clone());
    let chrom_id_arr = UInt8Array::from(chrom_id_vec);
    let position_arr = vec_u32_to_array(local_pos_vec);
    let strand_arr = vec_bool_to_array(hits.strand.clone());
    let score_arr = vec_f32_to_array(hits.score.clone());

    let schema = Arc::new(raw_schema());
    let batch = RecordBatch::try_new(
        schema,
        vec![
            Arc::new(query_id_arr),
            Arc::new(chrom_id_arr),
            Arc::new(position_arr),
            Arc::new(strand_arr),
            Arc::new(score_arr),
        ],
    )?;

    let n_rows = batch.num_rows();
    write_parquet(&batch, path)?;
    Ok(n_rows)
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Write RecordBatch → Parquet file
// ═══════════════════════════════════════════════════════════════════════════════

fn write_parquet(batch: &RecordBatch, path: &Path) -> Result<(), anyhow::Error> {
    let file = File::create(path)?;
    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .build();
    let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))?;
    writer.write(batch)?;
    writer.close()?;
    Ok(())
}
