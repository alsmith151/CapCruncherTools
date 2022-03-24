import logging
import pathlib
import subprocess
from . import cct


def fastq_deduplicate(
    infiles: tuple,
    output_prefix: str,
    output_compression: bool = False,
    n_cores: int = 1,
    **kwargs
):
    """Deduplicates paired FASTQ files.
    Input files in format: fq1_1,fq2_1,fq3_1 fq1_2,fq2_2,fq3_2
    """
    fq_in = {
        str(ii): [r1, r2]
        for ii, (r1, r2) in enumerate(
            zip(
                *[
                    files.split(",")
                    if "," in files
                    else [
                        files,
                    ]
                    for files in infiles
                ]
            )
        )
    }

    fq_out = {
        str(ii): [
            f"{output_prefix}{pathlib.Path(read).name.replace('.gz', '')}"
            for read in (r1, r2)
        ]
        for ii, (r1, r2) in enumerate(
            zip(
                *[
                    files.split(",")
                    if "," in files
                    else [
                        files,
                    ]
                    for files in infiles
                ]
            )
        )
    }

    #logging.debug(f"Input: {fq_in}")
    #logging.debug(f"Output: {fq_out}")

    deduplication_results = cct.deduplicate_fastq(
        fq_in,  # Fastq input list
        fq_out,  # Fastq output list
        True,  # Shuffle files
        output_compression,  # Compress output
    )

    #logging.debug(f"Deduplication results: {deduplication_results}")

    return deduplication_results
