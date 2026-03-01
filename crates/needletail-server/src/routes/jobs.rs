use std::sync::Arc;
use std::sync::atomic::Ordering;

use axum::extract::{Path, State};
use axum::http::StatusCode;
use axum::Json;
use serde_json::{json, Value};

use needletail_core::io::json::region_to_json;

use crate::jobs::JobStatus;
use crate::AppState;

/// Get job status and progress.
pub async fn status(
    State(state): State<Arc<AppState>>,
    Path(id): Path<String>,
) -> Result<Json<Value>, (StatusCode, Json<Value>)> {
    let job = state.jobs.get(&id).ok_or_else(|| {
        (
            StatusCode::NOT_FOUND,
            Json(json!({ "error": "job not found" })),
        )
    })?;

    let status = *job.status.lock().unwrap();
    let stage = job.progress.stage.lock().unwrap().clone();
    let current = job.progress.current.load(Ordering::Acquire);
    let total = job.progress.total.load(Ordering::Acquire);

    let status_str = match status {
        JobStatus::Pending => "pending",
        JobStatus::Running => "running",
        JobStatus::Completed => "completed",
        JobStatus::Failed => "failed",
        JobStatus::Cancelled => "cancelled",
    };

    let mut resp = json!({
        "id": job.id,
        "genome_id": job.genome_id,
        "status": status_str,
        "progress": {
            "stage": stage,
            "current": current,
            "total": total,
        }
    });

    if status == JobStatus::Failed {
        if let Some(ref err) = *job.error.lock().unwrap() {
            resp["error"] = json!(err);
        }
    }

    if status == JobStatus::Completed {
        if let Some(Ok(ref result)) = *job.result.lock().unwrap() {
            resp["summary"] = json!({
                "n_guides": result.guides.len(),
                "n_feature_tiles": result.feature_tiles.len(),
                "total_guides_scored": result.total_guides_scored,
            });
        }
    }

    Ok(Json(resp))
}

/// Stream job results as JSON array.
pub async fn stream(
    State(state): State<Arc<AppState>>,
    Path(id): Path<String>,
) -> Result<Json<Value>, (StatusCode, Json<Value>)> {
    let job = state.jobs.get(&id).ok_or_else(|| {
        (
            StatusCode::NOT_FOUND,
            Json(json!({ "error": "job not found" })),
        )
    })?;

    let status = *job.status.lock().unwrap();
    if status != JobStatus::Completed {
        return Err((
            StatusCode::BAD_REQUEST,
            Json(json!({ "error": format!("job status is {:?}, expected completed", status) })),
        ));
    }

    let result = job.result.lock().unwrap();
    match result.as_ref() {
        Some(Ok(lib_result)) => {
            let json_guides: Vec<Value> = lib_result.guides.iter().map(region_to_json).collect();
            Ok(Json(json!(json_guides)))
        }
        Some(Err(e)) => Err((
            StatusCode::INTERNAL_SERVER_ERROR,
            Json(json!({ "error": e })),
        )),
        None => Err((
            StatusCode::INTERNAL_SERVER_ERROR,
            Json(json!({ "error": "no result available" })),
        )),
    }
}

/// Cancel a running job.
pub async fn cancel(
    State(state): State<Arc<AppState>>,
    Path(id): Path<String>,
) -> Result<Json<Value>, (StatusCode, Json<Value>)> {
    if state.jobs.cancel(&id) {
        Ok(Json(json!({ "cancelled": true })))
    } else {
        Err((
            StatusCode::NOT_FOUND,
            Json(json!({ "error": "job not found" })),
        ))
    }
}
