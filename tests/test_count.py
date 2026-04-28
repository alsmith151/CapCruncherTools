from pathlib import Path

import pytest
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
