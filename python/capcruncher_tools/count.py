import pathlib

import pandas as pd
import polars as pl


def get_viewpoint(
    parquet: pathlib.Path,
    viewpoint: str,
    part: str | int | None = None,
    remove_exclusions: bool = False,
    remove_viewpoint: bool = False,
    subsample: float = 0,
    scan_low_memory: bool = False,
) -> pl.DataFrame:
    with pl.StringCache():

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
        counts = pd.DataFrame(counts.to_dicts())
    return counts


def count_viewpoint_pixels(
    parquet: str,
    viewpoint: str,
    remove_exclusions: bool = False,
    remove_viewpoint: bool = False,
    subsample: float = 0,
    low_memory: bool = False,
    partitions: list[str] | None = None,
) -> tuple[str, pd.DataFrame]:
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
        counts = counts.group_by(["bin1_id", "bin2_id"]).agg(pl.sum("count"))
        counts = pd.DataFrame(counts.to_dicts())

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
