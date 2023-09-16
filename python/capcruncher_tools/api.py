import os
import pathlib

import tempfile
import pandas as pd
import polars as pl
import tabulate
from loguru import logger as logging
from typing import Union, Tuple, List, Literal, Dict

from .capcruncher_tools import deduplicate, digest


def deduplicate_fastq(
    fastq1: List[str],
    fastq2: List[str],
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

    output_prefix_path = pathlib.Path(output_prefix)
    output_prefix_path.mkdir(parents=True, exist_ok=True)
    

    fastq_in = [(str(f1), str(f2)) for f1,f2 in zip(fastq1, fastq2)]
    fastq_out = [
        (
            str(output_prefix_path / pathlib.Path(f1).name),
            str(output_prefix_path / pathlib.Path(f2).name),
        )
        for f1, f2 in fastq_in
    ]

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
    fastqs: List[str] = None,
    output: str = "digested.fastq.gz",
    read_type: Literal["flashed", "pe"] = "pe",
    restriction_enzyme: str = "dpnii",
    minimum_slice_length: int = 18,
    sample_name: str = "sample",
) -> Dict[str, pl.DataFrame]:
    """
    Digest FASTQ files.

    Args:
        fastqs: List of FASTQ files.
        read_type: Read type.
        restriction_enzyme: Restriction enzyme.
        output: Output file name.
        minimum_slice_length: Minimum slice length.
        sample_name: Sample name.

    Returns:
        DataFrame with digestion stats.
    """
       

    (
        stats_read,
        stats_hist_unfilt,
        stats_hist_filt,
        stats_hist_length,
    ) = digest.digest_fastq(
        fastqs,
        output,
        restriction_enzyme,
        read_type,
        sample_name,
        minimum_slice_length,
    )

    return {
        "stats_read_level": stats_read,
        "stats_hist_unfilt": stats_hist_unfilt,
        "stats_hist_filt": stats_hist_filt,
        "stats_hist_length": stats_hist_length,
    }


def digest_genome(
    fasta: str,
    output: str = "digested.bed",
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
        fasta,
        restriction_enzyme,
        output,
        remove_recognition_site,
        minimum_slice_length,
        n_threads,
    )


def count_interactions(
    reporters: os.PathLike,
    output: os.PathLike = "CC_cooler.hdf5",
    remove_exclusions: bool = False,
    remove_viewpoint: bool = False,
    subsample: float = 0,
    fragment_map: os.PathLike = None,
    viewpoint_path: os.PathLike = None,
    n_cores: int = 1,
    assay: Literal["capture", "tri", "tiled"] = "capture",
    **kwargs,
) -> os.PathLike:
    """
    Counts interactions between the viewpoint and the rest of the genome.

    Args:
        reporters: Path to reporters file.
        output: Output file name.
        remove_exclusions: Remove excluded regions.
        remove_viewpoint: Remove capture regions.
        subsample: Subsample reads.
        fragment_map: Path to fragment map.
        viewpoint_path: Path to viewpoint file.
        n_cores: Number of cores.
        assay: Assay type.
        **kwargs: Additional arguments.
    Returns:
        Path to the generated cooler file.

    """

    import random
    import string

    import ray
    import tqdm
    import pyranges as pr

    import capcruncher.api.storage
    import capcruncher_tools.count

    # Extract viewpoint names and sizes from the supplied parquet file
    logging.info("Extracting viewpoint names and sizes")

    df = pd.read_parquet(reporters, engine="pyarrow", columns=["viewpoint"])
    viewpoints = df["viewpoint"].cat.categories.to_list()
    viewpoint_sizes = df["viewpoint"].value_counts()
    viewpoint_sizes_dict = viewpoint_sizes.to_dict()
    viewpoint_sizes_df = pd.DataFrame.from_dict(
        viewpoint_sizes_dict, orient="index", columns=["n_slices"]
    )
    viewpoint_sizes_df_tab = tabulate.tabulate(
        viewpoint_sizes_df, headers="keys", tablefmt="psql", showindex=True
    )

    logging.info(f"Number of viewpoints: {len(viewpoints)}")
    logging.info(f"Number of slices per viewpoint:\n {viewpoint_sizes_df_tab}")

    # Select a running mode
    if any([vp for vp in viewpoint_sizes_dict.values() if vp > 2e6]):
        logging.warning(
            "High number of slices per viewpoint detected. Switching to low memory mode"
        )
        low_memory = True

        # Extract the partitions used to generate the file
        df = pd.read_parquet(reporters, engine="pyarrow", columns=["bam"])
        partitions = df["bam"].cat.categories.to_list()

    else:
        low_memory = False

    # Start a Ray cluster
    ray.init(num_cpus=n_cores, ignore_reinit_error=True)

    # Store a reference to the bins table in the object store
    bins = pr.read_bed(fragment_map, as_df=True).rename(
        columns={
            "Chromosome": "chrom",
            "Start": "start",
            "End": "end",
            "Name": "name",
        }
    )
    bins_ref = ray.put(bins)

    # Fix reporters path
    reporters_path = pathlib.Path(reporters)
    if reporters_path.is_dir():
        reporters = str(reporters_path / "*.parquet")

    # Create a list of count futures
    futures = []
    for viewpoint in viewpoints:
        futures.append(
            capcruncher_tools.count.count_interactions.remote(
                parquet=f"{reporters}",
                viewpoint=viewpoint,
                remove_exclusions=remove_exclusions,
                remove_viewpoint=remove_viewpoint,
                subsample=subsample,
                low_memory=low_memory,
                partitions=partitions if low_memory else None,
            )
        )

    with tempfile.TemporaryDirectory() as tmpdir:
        cooler_refs = []
        for count_future in futures:
            cooler_refs.append(
                capcruncher_tools.count.make_cooler.remote(
                    output_prefix=str(
                        pathlib.Path(tmpdir)
                        / f"{''.join(random.choices(string.ascii_letters + string.digits, k=8))}.hdf5"
                    ),
                    future=count_future,
                    bins=bins_ref,
                    viewpoint_path=viewpoint_path,
                    assay=assay,
                )
            )

        with tqdm.tqdm(total=len(cooler_refs)) as pbar:
            coolers = []
            while cooler_refs:
                coolers_completed, cooler_refs = ray.wait(cooler_refs)
                for clr in coolers_completed:
                    clr = ray.get(clr)
                    coolers.append(clr.split("::")[0])
                    pbar.update(1)

        logging.info(f"Making final cooler at {output}")
        capcruncher.api.storage.merge_coolers(coolers, output=output)

    return output
