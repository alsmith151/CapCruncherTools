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

Dependency ranges are intentionally broad so `capcruncher-tools` can resolve in
CapCruncher pipeline environments without forcing the pipeline dependency set.
CapCruncher itself is expected to be supplied by the parent pipeline or CLI
environment.

## Current Tools

### FASTQ deduplication - fastq-deduplicate:

This tool takes paired FASTQ files and removes any duplicate fragments. 

#### Basic Usage:

```bash
capcruncher-tools fastq-deduplicate -1 <input1.fastq> -2 <input2.fastq> -o <output_prefix>
```


### Restriction digestion of FASTA - digest-genome:

This tool takes a FASTA file and a list of restriction enzymes and produces a list of fragments in BED format.

#### Basic Usage:

```bash
capcruncher-tools digest-genome -i <input.fasta> -o <output.bed> -r <recognition site> -p  <number of threads>
```


### Count restriction fragments - count:

This tool counts the number of interactions between a fragment (in silico digested read).

#### Basic Usage:

```bash
capcruncher-tools count <reporters> -f <fragments.bed> -v <viewpoints_path.bed> -o <output.hdf5> -p <number of threads>
```
