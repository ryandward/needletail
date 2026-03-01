//! JSON serialization for Region output.
//!
//! Converts Regions to JSON objects and supports streaming JSON array output.
//! NaN/Inf values are sanitized to null.

use std::io::Write;

use serde_json::{json, Map, Value};

use crate::models::region::{Region, TagValue};

/// Convert a Region to a JSON value.
pub fn region_to_json(r: &Region) -> Value {
    let mut obj = Map::new();

    obj.insert("chrom".into(), json!(r.chrom));
    obj.insert("start".into(), json!(r.start));
    obj.insert("end".into(), json!(r.end));
    obj.insert("strand".into(), json!(r.strand.as_str()));
    obj.insert("score".into(), sanitize_float(r.score));
    obj.insert("name".into(), json!(r.name));

    // Flatten tags into the object
    for (key, val) in &r.tags {
        obj.insert(key.clone(), tag_value_to_json(val));
    }

    Value::Object(obj)
}

/// Convert a TagValue to a JSON value.
fn tag_value_to_json(v: &TagValue) -> Value {
    match v {
        TagValue::Int(i) => json!(i),
        TagValue::Float(f) => sanitize_float(*f),
        TagValue::Str(s) => json!(s),
        TagValue::Bool(b) => json!(b),
    }
}

/// Sanitize a float: NaN and Inf become null.
fn sanitize_float(f: f64) -> Value {
    if f.is_finite() {
        json!(f)
    } else {
        Value::Null
    }
}

/// Write a JSON array of Regions to a writer, streaming one object at a time.
pub fn stream_json_array<W: Write>(regions: &[Region], writer: &mut W) -> std::io::Result<()> {
    writer.write_all(b"[")?;
    for (i, r) in regions.iter().enumerate() {
        if i > 0 {
            writer.write_all(b",")?;
        }
        let val = region_to_json(r);
        serde_json::to_writer(&mut *writer, &val)?;
    }
    writer.write_all(b"]")?;
    Ok(())
}

/// Serialize Regions to a JSON string.
pub fn regions_to_json_string(regions: &[Region]) -> String {
    let mut buf = Vec::new();
    stream_json_array(regions, &mut buf).expect("writing to Vec never fails");
    String::from_utf8(buf).expect("JSON is always valid UTF-8")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::region::Strand;

    #[test]
    fn test_region_to_json() {
        let r = Region::new("chr1", 100, 200)
            .with_strand(Strand::Forward)
            .with_name("guide1")
            .with_score(2.0)
            .with_tag("spacer", "ATCG")
            .with_tag("off_targets", 5i64);

        let j = region_to_json(&r);
        assert_eq!(j["chrom"], "chr1");
        assert_eq!(j["start"], 100);
        assert_eq!(j["end"], 200);
        assert_eq!(j["strand"], "+");
        assert_eq!(j["score"], 2.0);
        assert_eq!(j["spacer"], "ATCG");
        assert_eq!(j["off_targets"], 5);
    }

    #[test]
    fn test_nan_sanitization() {
        let r = Region::new("chr1", 0, 10).with_score(f64::NAN);
        let j = region_to_json(&r);
        assert!(j["score"].is_null());
    }

    #[test]
    fn test_stream_json_array() {
        let regions = vec![
            Region::new("chr1", 0, 10).with_name("a"),
            Region::new("chr1", 20, 30).with_name("b"),
        ];
        let json_str = regions_to_json_string(&regions);
        let parsed: Vec<Value> = serde_json::from_str(&json_str).unwrap();
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0]["name"], "a");
        assert_eq!(parsed[1]["name"], "b");
    }
}
