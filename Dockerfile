# syntax=docker/dockerfile:1
# ── Stage 1: Build ────────────────────────────────────────────────────────────
FROM rust:1.85-bookworm AS builder

WORKDIR /build

# Copy manifests first for dependency caching
COPY Cargo.toml Cargo.lock ./
COPY crates/needletail-core/Cargo.toml crates/needletail-core/Cargo.toml
COPY crates/needletail-py/Cargo.toml   crates/needletail-py/Cargo.toml
COPY crates/needletail-server/Cargo.toml crates/needletail-server/Cargo.toml

# Dummy sources so cargo can resolve the workspace and fetch deps
RUN mkdir -p crates/needletail-core/src && echo "" > crates/needletail-core/src/lib.rs \
 && mkdir -p crates/needletail-py/src   && echo "" > crates/needletail-py/src/lib.rs \
 && mkdir -p crates/needletail-server/src && echo "fn main() {}" > crates/needletail-server/src/main.rs

# Pre-build deps with BuildKit cache mounts (survives across builds)
RUN --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/build/target \
    cargo build --release -p needletail-server 2>/dev/null || true

# Copy real source and embedded assets
COPY crates/ crates/
COPY presets/ presets/

# Invalidate dummy build, then full release build (cached deps skip recompile)
RUN touch crates/needletail-core/src/lib.rs crates/needletail-server/src/main.rs
RUN --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/build/target \
    cargo build --release -p needletail-server \
 && cp target/release/needletail-server /usr/local/bin/needletail-server

# ── Stage 2: Runtime ─────────────────────────────────────────────────────────
FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/bin/needletail-server /usr/local/bin/

ENV BIND_ADDR=0.0.0.0:8002

EXPOSE 8002

CMD ["needletail-server"]
