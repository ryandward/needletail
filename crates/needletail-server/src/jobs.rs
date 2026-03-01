//! Background job management for pipeline execution.

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Arc;

use dashmap::DashMap;
use uuid::Uuid;

use needletail_core::models::preset::{CRISPRPreset, FeatureConfig};
use needletail_core::models::region::Region;
use needletail_core::pipeline::design::{self, LibraryResult, NullProgress, ProgressSink};

use crate::stores::genome_store::StoredGenome;

/// Job status.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JobStatus {
    Pending,
    Running,
    Completed,
    Failed,
    Cancelled,
}

/// Shared progress state for a running job.
pub struct JobProgress {
    pub stage: std::sync::Mutex<String>,
    pub current: AtomicUsize,
    pub total: AtomicUsize,
    pub cancelled: AtomicBool,
}

impl ProgressSink for JobProgress {
    fn report(&self, stage: &str, current: usize, total: usize) {
        *self.stage.lock().unwrap() = stage.to_string();
        self.current.store(current, Ordering::Release);
        self.total.store(total, Ordering::Release);
    }

    fn is_cancelled(&self) -> bool {
        self.cancelled.load(Ordering::Acquire)
    }
}

/// A completed or in-progress job.
pub struct Job {
    pub id: String,
    pub genome_id: String,
    pub status: std::sync::Mutex<JobStatus>,
    pub progress: Arc<JobProgress>,
    pub result: std::sync::Mutex<Option<Result<LibraryResult, String>>>,
    pub error: std::sync::Mutex<Option<String>>,
}

/// Thread-safe job manager.
pub struct JobManager {
    jobs: DashMap<String, Arc<Job>>,
}

impl JobManager {
    pub fn new() -> Self {
        JobManager {
            jobs: DashMap::new(),
        }
    }

    /// Submit a design_library job to run in the background.
    pub fn submit(
        &self,
        genome: Arc<StoredGenome>,
        preset: CRISPRPreset,
        feature_config: FeatureConfig,
    ) -> String {
        let job_id = Uuid::new_v4().to_string();
        let progress = Arc::new(JobProgress {
            stage: std::sync::Mutex::new("pending".into()),
            current: AtomicUsize::new(0),
            total: AtomicUsize::new(0),
            cancelled: AtomicBool::new(false),
        });

        let job = Arc::new(Job {
            id: job_id.clone(),
            genome_id: genome.id.clone(),
            status: std::sync::Mutex::new(JobStatus::Pending),
            progress: progress.clone(),
            result: std::sync::Mutex::new(None),
            error: std::sync::Mutex::new(None),
        });

        self.jobs.insert(job_id.clone(), job.clone());

        // Spawn blocking task
        tokio::task::spawn_blocking(move || {
            *job.status.lock().unwrap() = JobStatus::Running;

            let result = design::design_library(
                &genome.genome,
                &genome.index,
                genome.tier_small.as_ref(),
                genome.tier_large.as_ref(),
                &preset,
                &feature_config,
                &*progress,
            );

            match &result {
                Ok(_) => {
                    *job.status.lock().unwrap() = JobStatus::Completed;
                }
                Err(e) => {
                    if e == "cancelled" {
                        *job.status.lock().unwrap() = JobStatus::Cancelled;
                    } else {
                        *job.status.lock().unwrap() = JobStatus::Failed;
                        *job.error.lock().unwrap() = Some(e.clone());
                    }
                }
            }

            *job.result.lock().unwrap() = Some(result);
        });

        job_id
    }

    /// Get a job by ID.
    pub fn get(&self, id: &str) -> Option<Arc<Job>> {
        self.jobs.get(id).map(|v| v.clone())
    }

    /// Cancel a job.
    pub fn cancel(&self, id: &str) -> bool {
        if let Some(job) = self.jobs.get(id) {
            job.progress.cancelled.store(true, Ordering::Release);
            true
        } else {
            false
        }
    }

    /// List all job IDs with their statuses.
    pub fn list(&self) -> Vec<(String, String, JobStatus)> {
        self.jobs
            .iter()
            .map(|entry| {
                let job = entry.value();
                let status = *job.status.lock().unwrap();
                (job.id.clone(), job.genome_id.clone(), status)
            })
            .collect()
    }
}
