import pathlib

import pandas as pd
import polars as pl


def get_viewpoint(
    parquet: pathlib.Path,
    viewpoint: str,
    remove_exclusions: bool = False,
    remove_viewpoint: bool = False,
    subsample: float = 0,
    low_memory: bool = False,
):
    df = pl.scan_parquet(parquet, low_memory=low_memory).filter(
        pl.col("viewpoint") == viewpoint
    )

    if remove_viewpoint:
        df = df.filter(pl.col("capture_count") == 0)

    if remove_exclusions:
        df = df.filter(pl.col("exclusion") != pl.col("viewpoint"))

    df = df.select(["parent_id", "restriction_fragment"])

    return df.collect()


def get_counts(df: pl.DataFrame):
    from capcruncher_tools import interactions

    counts = interactions.count_interactions(df)
    return counts.to_pandas()


def count_viewpoint_pixels(
    parquet: pathlib.Path,
    viewpoint: str,
    remove_exclusions: bool = False,
    remove_viewpoint: bool = False,
    subsample: float = 0,
    low_memory: bool = False,
    partitions: list[str] | None = None,
):
    """Count pixels for one CapCruncher viewpoint."""
    if not low_memory or not partitions:
        df = get_viewpoint(
            parquet,
            viewpoint,
            remove_exclusions=remove_exclusions,
            remove_viewpoint=remove_viewpoint,
            subsample=subsample,
            low_memory=low_memory,
        )
        return viewpoint, get_counts(df)

    counts = []
    for partition in partitions:
        df = (
            pl.scan_parquet(parquet, low_memory=True)
            .filter((pl.col("viewpoint") == viewpoint) & (pl.col("bam") == partition))
        )

        if remove_viewpoint:
            df = df.filter(pl.col("capture_count") == 0)

        if remove_exclusions:
            df = df.filter(pl.col("exclusion") != pl.col("viewpoint"))

        df = df.select(["parent_id", "restriction_fragment"]).collect()
        if not df.is_empty():
            counts.append(get_counts(df))

    if not counts:
        return viewpoint, get_counts(
            pl.DataFrame(
                schema={
                    "parent_id": pl.Int64,
                    "restriction_fragment": pl.Int64,
                }
            )
        )

    return (
        viewpoint,
        pd.concat(counts).groupby(["bin1_id", "bin2_id"], as_index=False).sum(),
    )
