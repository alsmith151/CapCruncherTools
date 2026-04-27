# CapCruncherTools

## Overview

A collection of Rust tools to speed up the Python-based functionality of the CapCruncher project. Python bindings have been generated to allow for easy integration into the existing Python code base.

## Development

This project uses [uv](https://docs.astral.sh/uv/) for Python dependency management and environment setup. Dependencies are declared in `pyproject.toml` and resolved in `uv.lock`.

Create or update the local environment:

```bash
uv sync --dev
```

Build the Rust extension into the uv environment:

```bash
uv run maturin develop --release
```

Run the test suite:

```bash
uv run pytest
```

Refresh locked dependencies after editing `pyproject.toml`:

```bash
uv lock --upgrade
```

Runtime dependency ranges are intentionally broad so `capcruncher-tools` can
install inside CapCruncher pipeline environments without forcing shared
libraries such as Click, pandas, or Polars to a single version. Build and
development tools remain pinned for reproducible local builds. CapCruncher
itself is expected to be supplied by the parent pipeline or CLI environment.

## Current Accelerators

### FASTQ deduplication

This tool takes paired FASTQ files and removes any duplicate fragments. 

Use it through the CapCruncher CLI or Python API; this package provides the
compiled acceleration layer.


### Restriction digestion of FASTA

This tool takes a FASTA file and a list of restriction enzymes and produces a list of fragments in BED format.

Use it through the CapCruncher CLI or Python API; this package does not install a
separate `capcruncher-tools` command.


### Count restriction fragments

This tool counts the number of interactions between a fragment (in silico digested read).

Use it through CapCruncher's interaction-counting API; cooler creation and
workflow orchestration live in CapCruncher.
