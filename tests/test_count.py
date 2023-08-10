import os
import tempfile
import shutil
import pandas as pd
import pyranges as pr
from capcruncher_tools.count import count_interactions, make_cooler
from capcruncher_tools.cli import cli
import click.testing
import ray
import pytest

@pytest.fixture(scope="module")
def data_path():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "count")
    return data_dir

@pytest.mark.parametrize("assay", ["capture", "tri", "tiled"])
def test_count(data_path, assay):
    # Set up test data
    reporters = os.path.join(data_path, "SAMPLE-A_REP2.parquet")
    fragment_map = os.path.join(data_path, "bins.bed.gz")
    viewpoint_path = os.path.join(data_path, "viewpoints.bed")
    output = os.path.join(data_path, "counts.hdf5")

    # Set up temporary directory for output
    with tempfile.TemporaryDirectory() as tmpdir:
        # Run count function
        result = click.testing.CliRunner().invoke(
            cli,
            [
                "count",
                reporters,
                "-o",
                os.path.join(tmpdir, "counts.hdf5"),
                "-f",
                fragment_map,
                "-v",
                viewpoint_path,
                "--assay",
                assay,
            ],
        )

        if result.exit_code != 0:
            print(result.stdout)
        assert result.exit_code == 0

        # Check that output file exists
        assert os.path.exists(os.path.join(tmpdir, "counts.hdf5"))