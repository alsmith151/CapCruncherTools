import pathlib
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
