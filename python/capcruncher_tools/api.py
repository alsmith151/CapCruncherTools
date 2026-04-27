from __future__ import annotations

import os
from pathlib import Path
from typing import Sequence

import pandas as pd

from .capcruncher_tools import deduplicate, digest


def deduplicate_fastq(
    fastq1: Sequence[str | os.PathLike[str]],
    fastq2: Sequence[str | os.PathLike[str]],
    output_prefix: str | os.PathLike[str],
    sample_name: str = "sampleX",
    shuffle: bool = False,
) -> pd.DataFrame:
    """Deduplicate paired FASTQ files and return CapCruncher-style statistics."""
    fq_input = [(os.fspath(fq1), os.fspath(fq2)) for fq1, fq2 in zip(fastq1, fastq2)]
    fq_output = [
        (
            os.fspath(output_prefix) + os.path.basename(fq1),
            os.fspath(output_prefix) + os.path.basename(fq2),
        )
        for fq1, fq2 in fq_input
    ]

    if fq_output:
        Path(fq_output[0][0]).parent.mkdir(parents=True, exist_ok=True)

    deduplication_results = deduplicate.fastq_deduplicate(
        fq_input,
        fq_output,
        shuffle,
    )

    return (
        pd.Series(deduplication_results)
        .to_frame("stat")
        .reset_index()
        .rename(columns={"index": "stat_type"})
        .assign(
            read_number=0,
            read_type="pe",
            stage="deduplication",
            sample=sample_name,
        )
    )


def digest_genome(
    fasta: str | os.PathLike[str],
    output: str | os.PathLike[str],
    restriction_enzyme: str,
    remove_recognition_site: bool = True,
    minimum_slice_length: int = 18,
    n_threads: int = 1,
) -> None:
    """Digest a FASTA file to BED records for CapCruncher."""
    digest.digest_fasta(
        os.fspath(fasta),
        restriction_enzyme,
        os.fspath(output),
        remove_recognition_site,
        minimum_slice_length,
        n_threads,
    )
