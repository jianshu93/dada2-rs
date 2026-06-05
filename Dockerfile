# ---- Build stage ----
FROM rust:1-bookworm AS builder

# Optional pre-computed version string (e.g. "0.1.0-deadbeef" or "0.1.0").
# When unset, build.rs falls back to the bare Cargo.toml version because
# .git is not part of the Docker context.
ARG DADA2_RS_VERSION_FULL
ENV DADA2_RS_VERSION_FULL=${DADA2_RS_VERSION_FULL}

# Build tuning knobs. Defaults reproduce the portable build:
#   CARGO_PROFILE=release         -> target/release/dada2-rs
#   RUSTFLAGS=""                  -> baseline x86-64, runs anywhere
# The Docker workflow overrides these to produce an optimized image variant,
# e.g. CARGO_PROFILE=release-native + RUSTFLAGS="-C target-cpu=x86-64-v3".
# Note: a custom profile name (release-native) emits to target/<profile>/,
# while the built-in "release" emits to target/release/ — target/${PROFILE}/
# resolves correctly for both.
ARG CARGO_PROFILE=release
ARG RUSTFLAGS=""
ENV RUSTFLAGS=${RUSTFLAGS}

WORKDIR /build
COPY Cargo.toml Cargo.lock build.rs ./
COPY src/ src/

RUN cargo build --profile "${CARGO_PROFILE}" \
    && strip "target/${CARGO_PROFILE}/dada2-rs" \
    && cp "target/${CARGO_PROFILE}/dada2-rs" /tmp/dada2-rs

# ---- Runtime stage ----
FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /tmp/dada2-rs /usr/local/bin/dada2-rs

ENTRYPOINT ["dada2-rs"]
