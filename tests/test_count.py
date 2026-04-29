from pathlib import Path

import polars as pl
import pytest
from capcruncher_tools import interactions
from capcruncher_tools.count import count_viewpoint_pixels


@pytest.fixture(scope="module")
def data_path():
    return Path(__file__).resolve().parent / "count"


def test_count(data_path):
    reporters = data_path / "SAMPLE-A_REP2.parquet"
    viewpoint, counts = count_viewpoint_pixels(
        parquet=str(reporters / "*.parquet"),
        viewpoint="Slc25A37",
    )

    assert viewpoint == "Slc25A37"
    assert {"bin1_id", "bin2_id", "count"} <= set(counts.columns)
    assert not counts.empty


def test_count_interactions_schema_contract_uses_unsigned_parent_ids():
    df = pl.DataFrame(
        {
            "parent_id": [1, 1],
            "restriction_fragment": [169686, 169744],
        },
        schema={
            "parent_id": pl.UInt64,
            "restriction_fragment": pl.Int64,
        },
    )

    counts = interactions.count_interactions(df)

    assert df.schema["parent_id"] == pl.UInt64
    assert counts.schema == {
        "bin1_id": pl.Int64,
        "bin2_id": pl.Int64,
        "count": pl.Int32,
    }


def test_count_interactions_preserves_unsigned_parent_ids_above_i64():
    df = pl.DataFrame(
        {
            "parent_id": [7, 7, 2**63 + 5, 2**63 + 5],
            "restriction_fragment": [169686, 169744, 169686, 169744],
        },
        schema={
            "parent_id": pl.UInt64,
            "restriction_fragment": pl.Int64,
        },
    )

    counts = interactions.count_interactions(df).sort(["bin1_id", "bin2_id"])

    assert counts.to_dicts() == [
        {"bin1_id": 169686, "bin2_id": 169744, "count": 2},
    ]
