use std::sync::Arc;

use axum::extract::State;
use axum::http::StatusCode;
use axum::Json;
use serde::Deserialize;
use serde_json::{json, Value};

use needletail_core::models::preset::{CRISPRPreset, FeatureConfig};

use crate::AppState;

/// List available methods.
pub async fn catalog() -> Json<Value> {
    Json(json!([
        {
            "name": "design_library",
            "description": "Design a CRISPR guide library for a genome",
            "parameters": {
                "genome_id": "string (required)",
                "preset": "string (optional, default: 'spcas9')",
                "feature_config": "string (optional, default: 'saccer3')",
                "pam": "string (optional, overrides preset PAM)",
                "spacer_len": "integer (optional, overrides preset spacer_len)",
                "mismatches": "integer (optional, default: 0)"
            }
        }
    ]))
}

#[derive(Deserialize)]
pub struct DesignRequest {
    pub genome_id: String,
    pub preset: Option<String>,
    pub feature_config: Option<String>,
    pub pam: Option<String>,
    pub spacer_len: Option<usize>,
    pub mismatches: Option<u8>,
}

/// Dispatch a design_library job.
pub async fn design_library(
    State(state): State<Arc<AppState>>,
    Json(req): Json<DesignRequest>,
) -> Result<Json<Value>, (StatusCode, Json<Value>)> {
    let genome = state.genomes.get(&req.genome_id).ok_or_else(|| {
        (
            StatusCode::NOT_FOUND,
            Json(json!({ "error": "genome not found" })),
        )
    })?;

    // Build preset
    let mut preset = match req.preset.as_deref() {
        Some("cas12a") => CRISPRPreset::cas12a(),
        _ => CRISPRPreset::spcas9(),
    };
    if let Some(pam) = req.pam {
        preset.pam = pam;
    }
    if let Some(sl) = req.spacer_len {
        preset.spacer_len = sl;
    }
    if let Some(mm) = req.mismatches {
        preset.mismatches = mm;
    }

    // Build feature config
    let feature_config = FeatureConfig::saccer3();

    let job_id = state.jobs.submit(genome, preset, feature_config);

    Ok(Json(json!({ "job_id": job_id })))
}
