import pathlib
from typing import List, Literal, Tuple, Union

import capcruncher.api as cc
import pandas as pd
import polars as pl
import ray

pl.enable_string_cache(True)


def get_viewpoint(
    parquet: pathlib.Path,
    viewpoint: str,
    part: Union[str, int] = None,
    remove_exclusions: bool = False,
    remove_viewpoint: bool = False,
    subsample: float = 0,
    scan_low_memory: bool = False,
) -> pl.DataFrame:
    if not part:
        df = pl.scan_parquet(parquet, low_memory=scan_low_memory).filter(
            pl.col("viewpoint") == viewpoint
        )

    else:
        df = pl.scan_parquet(parquet, low_memory=scan_low_memory).filter(
            (pl.col("viewpoint") == viewpoint) & (pl.col("bam") == part)
        )

    if remove_viewpoint:
        df = df.filter(pl.col("capture_count") == 0)

    if remove_exclusions:
        df = df.filter(pl.col("exclusion") != pl.col("viewpoint"))

    df = df.select(["parent_id", "restriction_fragment"])

    return df.collect()


def get_counts(df: pl.DataFrame, as_pandas: bool = True) -> pd.DataFrame:
    from capcruncher_tools import interactions

    counts = interactions.count_interactions(df)
    if as_pandas:
        counts = counts.to_pandas()
    return counts


@ray.remote
def count_interactions(
    parquet: str,
    viewpoint: str,
    remove_exclusions: bool = False,
    remove_viewpoint: bool = False,
    subsample: float = 0,
    low_memory: bool = False,
    partitions: List[str] = None,
) -> Tuple[str, pd.DataFrame]:
    from .count import get_counts, get_viewpoint

    if low_memory:
        # Check that partitions are specified
        assert (
            partitions is not None
        ), "partitions must be specified when using low_memory"

        # Get counts for each partition
        counts = []
        for partition in partitions:
            df = get_viewpoint(
                parquet=parquet,
                viewpoint=viewpoint,
                remove_exclusions=remove_exclusions,
                remove_viewpoint=remove_viewpoint,
                subsample=subsample,
                scan_low_memory=True,
                part=partition,
            )

            count = get_counts(df, as_pandas=False)
            counts.append(count)

        # Combine counts
        counts = pl.concat(counts)
        counts = counts.groupby(["bin1_id", "bin2_id"]).agg(pl.sum("count"))
        counts = counts.to_pandas()

    else:
        df = get_viewpoint(
            parquet=parquet,
            viewpoint=viewpoint,
            remove_exclusions=remove_exclusions,
            remove_viewpoint=remove_viewpoint,
            subsample=subsample,
            scan_low_memory=False,
        )
        counts = get_counts(df)

    return (viewpoint, counts)


@ray.remote
def make_cooler(
    output_prefix: str,
    future: "ray.ObjectRef",
    bins: pd.DataFrame,
    viewpoint_path: str,
    assay: Literal["capture", "tri", "tiled"],
    **kwargs,
) -> str:
    viewpoint_name, counts = future

    return cc.storage.create_cooler_cc(
        output_prefix=output_prefix,
        pixels=counts,
        bins=bins,
        viewpoint_name=viewpoint_name,
        viewpoint_path=viewpoint_path,
        assay=assay,
        **kwargs,
    )
