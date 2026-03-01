//! Needletail HTTP Server — standalone Rust binary serving the CRISPR design API.
//!
//! Replaces SeqChain's Python HTTP layer with Axum.

mod jobs;
mod routes;
mod stores;

use std::sync::Arc;

use axum::routing::{delete, get, post};
use axum::Router;
use tower_http::cors::CorsLayer;

use crate::jobs::JobManager;
use crate::stores::genome_store::GenomeStore;

/// Shared application state.
pub struct AppState {
    pub genomes: GenomeStore,
    pub jobs: JobManager,
}

#[tokio::main]
async fn main() {
    let state = Arc::new(AppState {
        genomes: GenomeStore::new(),
        jobs: JobManager::new(),
    });

    let app = Router::new()
        // Health
        .route("/api/health", get(routes::health::health))
        // Genomes
        .route("/api/genomes/upload", post(routes::genomes::upload))
        .route("/api/genomes/", get(routes::genomes::list))
        .route("/api/genomes/{id}", get(routes::genomes::get))
        // Methods
        .route("/api/methods", get(routes::methods::catalog))
        .route(
            "/api/methods/design_library",
            post(routes::methods::design_library),
        )
        // Jobs
        .route("/api/jobs/{id}", get(routes::jobs::status))
        .route("/api/jobs/{id}/stream", get(routes::jobs::stream))
        .route("/api/jobs/{id}", delete(routes::jobs::cancel))
        .layer(CorsLayer::permissive())
        .with_state(state);

    let bind = "0.0.0.0:3000";
    println!("needletail-server listening on {}", bind);

    let listener = tokio::net::TcpListener::bind(bind).await.unwrap();
    axum::serve(listener, app).await.unwrap();
}
