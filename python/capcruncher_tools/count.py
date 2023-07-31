import pathlib
import polars as pl
import ray
from typing import Tuple
import pandas as pd
from capcruncher.api import storage


def get_viewpoint(
    parquet: pathlib.Path,
    viewpoint: str,
    remove_exclusions: bool = False,
    remove_viewpoint: bool = False,
    subsample: float = 0,
    low_memory: bool = False,
) -> pl.DataFrame:
    df = pl.scan_parquet(parquet, low_memory=low_memory).filter(
        pl.col("viewpoint") == viewpoint
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
) -> Tuple[str, pd.DataFrame]:
    from .count import get_viewpoint, get_counts

    df = get_viewpoint(
        parquet,
        viewpoint,
        remove_exclusions,
        remove_viewpoint,
        subsample,
        low_memory=low_memory,
    )

    if low_memory:
        # Break up the dataframe into chunks of 1e6 rows
        # and aggregate the counts for each chunk
        counts = []

        # Iterate over 1e6 row chunks
        for i in range(int(df.shape[0] / 1e6) + 1):
            start = int(i * 1e6)
            end = int(min(((i + 1) * 1e6), df.shape[0]))
            # Get counts for chunk
            counts.append(get_counts(df[start:end], as_pandas=False))

        counts = pl.concat(counts)
        counts = counts.groupby(["bin1_id", "bin2_id"]).agg(pl.sum("count"))
        counts = counts.to_pandas()

    else:
        counts = get_counts(df)

    return (viewpoint, counts)


@ray.remote
def make_cooler(
    output_prefix: str,
    counts: pd.DataFrame,
    bins: pd.DataFrame,
    viewpoint_name: str,
    viewpoint_path: str,
    **kwargs,
) -> str:
    return storage.create_cooler_cc(
        output_prefix=output_prefix,
        pixels=counts,
        bins=bins,
        viewpoint_name=viewpoint_name,
        viewpoint_path=viewpoint_path,
        **kwargs,
    )
