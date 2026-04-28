from __future__ import annotations

from os import PathLike
from pathlib import Path
from typing import Literal, Sequence

import pandas as pd
from loguru import logger as logging

from .capcruncher_tools import deduplicate, digest


class DigestFastqStats:
    """Small wrapper for Rust digestion stats returned to CapCruncher."""

    def __init__(self, data):
        self.data = data

    def model_dump_json(self) -> str:
        import json

        return json.dumps(self.data)


def deduplicate_fastq(
    fastq1: Sequence[str | PathLike[str]],
    fastq2: Sequence[str | PathLike[str]],
    output_prefix: str = "deduplicated_",
    sample_name: str = "sample",
    shuffle: bool = False,
) -> pd.DataFrame:
    """
    Deduplicate FASTQ files.

    Args:
        fastq1: List of FASTQ files (R1).
        fastq2: List of FASTQ files (R2).
        output: Output file name.
        sample_name: Sample name.
        shuffle: Shuffle reads before deduplication.

    Returns:
        DataFrame with deduplicated read stats.

    """

    if len(fastq1) != len(fastq2):
        raise ValueError("Number of FASTQ files in R1 and R2 must be equal")

    fastq_in = [
        (str(Path(f1)), str(Path(f2))) for f1, f2 in zip(fastq1, fastq2)
    ]

    # Create output file names
    fastq_out = []
    for fqs in fastq_in:
        fq_pair = []
        for fq in fqs:
            fq_name = Path(fq).name
            fq_pair.append(f"{output_prefix}{fq_name}")
        fastq_out.append(tuple(fq_pair))

    if fastq_out:
        Path(fastq_out[0][0]).parent.mkdir(parents=True, exist_ok=True)

    logging.info("Deduplicating FASTQ files")
    deduplication_results = deduplicate.fastq_deduplicate(fastq_in, fastq_out, shuffle)

    logging.info("Preparing deduplication stats")
    df_stats = (
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

    return df_stats


def digest_fastq(
    fastqs: Sequence[str | PathLike[str]] | None = None,
    output: str = "digested.fastq.gz",
    read_type: Literal["flashed", "pe"] = "pe",
    restriction_site: str = "dpnii",
    minimum_slice_length: int = 18,
    sample_name: str = "sample",
):
    """
    Digest FASTQ files.

    Args:
        fastqs: List of FASTQ files.
        read_type: Read type.
        restriction_site: Restriction enzyme site.
        output: Output file name.
        minimum_slice_length: Minimum slice length.
        sample_name: Sample name.

    Returns:
        DataFrame with digestion stats.
    """
    stats = digest.digest_fastq(
        [str(Path(fastq)) for fastq in fastqs or []],
        str(Path(output)),
        restriction_site,
        read_type.capitalize(),
        sample_name,
        minimum_slice_length,
    )
    return DigestFastqStats(stats)


def digest_genome(
    fasta: str | PathLike[str],
    output: str | PathLike[str] = "digested.bed",
    restriction_enzyme: str = "DpnII",
    remove_recognition_site: bool = True,
    minimum_slice_length: int = 18,
    n_threads: int = 1,
):
    """
    Digest genome.

    Args:
        fasta: FASTA file.
        output: Output file name.
        restriction_enzyme: Restriction enzyme.

    Returns:
        DataFrame with digestion stats.
    """

    logging.info("Digesting genome")
    digest.digest_fasta(
        str(Path(fasta)),
        restriction_enzyme,
        str(Path(output)),
        remove_recognition_site,
        minimum_slice_length,
        n_threads,
    )
