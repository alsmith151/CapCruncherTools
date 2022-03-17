
import logging
import pathlib
from . import cct


def fastq_deduplicate(infiles: tuple, output_prefix: str):
    """Deduplicates paired FASTQ files.
    Input files in format: fq1_1,fq2_1,fq3_1 fq1_2,fq2_2,fq3_2
    """
    fq_in = {
        str(ii): [r1, r2]
        for ii, (r1, r2) in enumerate(
            zip(*[files.split(",") if "," in files else [files, ] for files in infiles])
        )
    }

    fq_out = {
        str(ii): [f"{output_prefix}{pathlib.Path(read).name}" for read in (r1, r2)]
        for ii, (r1, r2) in enumerate(
            zip(*[files.split(",") if "," in files else [files,] for files in infiles])
        )
    }

    logging.debug(f"Input: {fq_in}")
    logging.debug(f"Output: {fq_out}")

    deduplication_results = cct.deduplicate_fastq(fq_in, fq_out, True)
    
    logging.debug(f"Deduplication results: {deduplication_results}")

    return deduplication_results