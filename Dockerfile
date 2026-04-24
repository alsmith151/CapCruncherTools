# syntax=docker/dockerfile:1.7

ARG PYTHON_VERSION=3.12

FROM python:${PYTHON_VERSION}-slim-bookworm AS builder

ENV CARGO_HOME=/usr/local/cargo \
    PATH=/usr/local/cargo/bin:/root/.local/bin:$PATH \
    UV_CACHE_DIR=/tmp/uv-cache \
    UV_LINK_MODE=copy

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        pkg-config \
    && rm -rf /var/lib/apt/lists/* \
    && curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
        | sh -s -- -y --profile minimal --default-toolchain stable \
    && curl -LsSf https://astral.sh/uv/install.sh | sh

WORKDIR /app

COPY pyproject.toml uv.lock Cargo.toml ./
COPY src ./src
COPY python ./python
COPY README.md LICENSE ./

RUN --mount=type=cache,target=/tmp/uv-cache \
    --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/usr/local/cargo/git \
    uv build --wheel --out-dir /dist

RUN --mount=type=cache,target=/tmp/uv-cache \
    uv export --locked --no-dev --no-emit-project --no-hashes --output-file /tmp/requirements.txt \
    && uv pip install --prefix=/install --strict --requirements /tmp/requirements.txt \
    && uv pip install --prefix=/install --no-deps /dist/*.whl

FROM python:${PYTHON_VERSION}-slim-bookworm AS runtime

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /install /usr/local

ENTRYPOINT ["capcruncher-tools"]
CMD ["--help"]
