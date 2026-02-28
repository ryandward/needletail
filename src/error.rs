/// Errors raised by the FM-Index layer.
#[derive(Debug)]
pub enum SearchError {
    Io(std::io::Error),
    Other(anyhow::Error),
}

impl std::fmt::Display for SearchError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SearchError::Io(e) => write!(f, "IO error: {}", e),
            SearchError::Other(e) => write!(f, "Search error: {}", e),
        }
    }
}

impl std::error::Error for SearchError {}

impl From<std::io::Error> for SearchError {
    fn from(e: std::io::Error) -> Self {
        SearchError::Io(e)
    }
}

impl From<anyhow::Error> for SearchError {
    fn from(e: anyhow::Error) -> Self {
        SearchError::Other(e)
    }
}
