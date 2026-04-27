import os
import pytest
from capcruncher_tools.count import count_viewpoint_pixels

@pytest.fixture(scope="module")
def data_path():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "count")
    return data_dir

def test_count(data_path):
    reporters = os.path.join(data_path, "SAMPLE-A_REP2.parquet")
    viewpoint, counts = count_viewpoint_pixels(
        parquet=os.path.join(reporters, "*.parquet"),
        viewpoint="Slc25A37",
    )

    assert viewpoint == "Slc25A37"
    assert {"bin1_id", "bin2_id", "count"} <= set(counts.columns)
    assert not counts.empty
