import logging
import pathlib
from typing import List, Tuple

import capcruncher_tools

def fastq_deduplicate(
    infiles: Tuple[Tuple[str, str]],
    outfiles: Tuple[Tuple[str, str]],
    error_rate: float = 0.01,
    shuffle: bool = False,
):

    stats = capcruncher_tools.fastq_deduplicate(
        infiles,
        outfiles,
        shuffle,
        error_rate)
    
    return stats