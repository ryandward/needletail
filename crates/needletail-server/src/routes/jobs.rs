//! Job polling and streaming endpoints.
//!
//! ```text
//! GET    /api/jobs/:id        → status + progress
//! GET    /api/jobs/:id/stream → result stream (JSON)
//! DELETE /api/jobs/:id        → cancel (best-effort); always 200
//! ```

use std::sync::Arc;
use std::sync::atomic::Ordering;

use axum::body::Body;
use axum::extract::{Path, State};
use axum::http::{HeaderValue, StatusCode, header};
use axum::response::Response;
use axum::Json;
use serde_json::{json, Value};
use tokio_util::io::ReaderStream;

use crate::jobs::JobStatus;
use crate::AppState;

/// Poll job status and progress.
pub async fn status(
    State(state): State<Arc<AppState>>,
    Path(id): Path<String>,
) -> Result<Json<Value>, (StatusCode, Json<Value>)> {
    let job = state.jobs.get(&id).ok_or_else(|| {
        (
            StatusCode::NOT_FOUND,
            Json(json!({ "error": "Job not found" })),
        )
    })?;

    let status = *job.status.lock().unwrap();
    let elapsed = job.elapsed_secs();

    // Read step/stage/items from progress
    let step = job.progress.step.lock().unwrap().clone();
    let stage = job.progress.step_stage.lock().unwrap().clone();
    let items_complete = job.progress.items_complete.load(Ordering::Acquire);
    let items_total = job.progress.items_total.load(Ordering::Acquire);

    // Items object: only present when items_total > 0
    let items = if items_total > 0 {
        json!({ "complete": items_complete, "total": items_total })
    } else {
        Value::Null
    };

    // pct_complete: per-step, derived from items when available
    let pct: Option<f64> = if items_total > 0 {
        Some((items_complete as f64 / items_total as f64 * 1000.0).round() / 1000.0)
    } else {
        None
    };

    // rate_per_sec: live throughput derived from the items counter.
    // This is the actual processing rate of whichever step is currently
    // running set_items() — scoring chunks or annotation drain.
    let rate: Option<f64> = match (elapsed, items_complete) {
        (Some(e), ic) if e > 0.0 && ic > 0 => Some((ic as f64 / e * 10.0).round() / 10.0),
        _ => None,
    };

    let mut resp = json!({
        "status": status.as_str(),
        "step": step,
        "stage": stage,
        "items": items,
        "progress": {
            "pct_complete": pct,
            "rate_per_sec": rate,
        },
        "error": Value::Null,
    });

    if status == JobStatus::Failed {
        if let Some(ref err) = *job.error.lock().unwrap() {
            resp["error"] = json!(err);
        }
    }

    Ok(Json(resp))
}

/// Stream the completed job result (one-shot).
///
/// Called once by GenomeHub after status reaches "complete".
/// Streams the Parquet file directly from disk via an open fd.
///
/// **Lifecycle**: open fd → unlink path → stream from orphaned inode →
/// fd drops → OS reclaims disk blocks.  No timers, no heuristics.
///
/// This endpoint is **consumed on read**: the first successful call
/// takes ownership of the result.  Subsequent calls return 410 Gone.
pub async fn stream(
    State(state): State<Arc<AppState>>,
    Path(id): Path<String>,
) -> Result<Response, (StatusCode, Json<Value>)> {
    let job = state.jobs.get(&id).ok_or_else(|| {
        (
            StatusCode::NOT_FOUND,
            Json(json!({ "error": "Job not found" })),
        )
    })?;

    let status = *job.status.lock().unwrap();
    if status != JobStatus::Complete {
        return Err((
            StatusCode::ACCEPTED,
            Json(json!({ "error": format!("Job not complete (status: {})", status.as_str()) })),
        ));
    }

    // Take the result path — one-shot: first caller wins.
    let result_path = job.result_path.lock().unwrap().take();
    let taken_result = job.result.lock().unwrap().take();

    match (result_path, taken_result) {
        (Some(path), Some(Ok(_lib_result))) => {
            // Open the file descriptor FIRST — we need the fd before unlinking.
            let file = tokio::fs::File::open(&path).await.map_err(|e| {
                (
                    StatusCode::INTERNAL_SERVER_ERROR,
                    Json(json!({ "error": format!("Failed to open result file: {}", e) })),
                )
            })?;

            // Unlink the directory entry immediately.  The VFS keeps the
            // inode alive as long as our open fd exists.  When the
            // ReaderStream drops (response complete or client disconnect),
            // the fd closes and the OS reclaims the disk blocks.
            // No timers.  No races.  Pure VFS reference counting.
            tokio::fs::remove_file(&path).await.ok();

            let stream = ReaderStream::new(file);
            let body = Body::from_stream(stream);

            // Evict the job from the store.
            drop(job);
            state.jobs.remove(&id);

            let mut response = Response::new(body);
            response.headers_mut().insert(
                header::CONTENT_TYPE,
                HeaderValue::from_static("application/vnd.apache.parquet"),
            );
            Ok(response)
        }
        (_, Some(Err(e))) => Err((
            StatusCode::INTERNAL_SERVER_ERROR,
            Json(json!({ "error": e })),
        )),
        // Result already consumed by a prior /stream call, or
        // never materialized (cancelled before completion).
        _ => Err((
            StatusCode::GONE,
            Json(json!({ "error": "Result already consumed or never produced. This endpoint is one-shot: the stream is consumed upon first read." })),
        )),
    }
}

/// Cancel a running job.
///
/// Best-effort: if the job has already completed or never existed,
/// returns 200 anyway.
pub async fn cancel(
    State(state): State<Arc<AppState>>,
    Path(id): Path<String>,
) -> Json<Value> {
    state.jobs.cancel(&id);
    Json(json!({ "ok": true }))
}
