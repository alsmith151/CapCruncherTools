import pytest
import os
import click.testing
import pandas as pd

from capcruncher_tools.cli import cli


@pytest.fixture(scope="module")
def data_path():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "fastq_deduplicate")
    return data_dir


@pytest.mark.parametrize(
    "infiles,prefix,n_duplicates_expected",
    [
        (
            ("duplicated_1.fastq.gz", "duplicated_2.fastq.gz"),
            "out_",
            538,
        )
    ],
)
def test_fastq_duplicate_removal(
    data_path, tmpdir, infiles, prefix, n_duplicates_expected
):

    infiles_paths = [os.path.join(data_path, fn) for fn in infiles]
    out_prefix = os.path.join(tmpdir, prefix)

    result = click.testing.CliRunner().invoke(
        cli,
        [
            "fastq-deduplicate",
            "-1",
            infiles_paths[0],
            "-2",
            infiles_paths[1],
            "-o",
            out_prefix,
        ],
    )
    assert result.exit_code == 0

    # Check that the output files exist
    outfiles = [
        out_prefix + os.path.basename(infiles_paths[0]),
        out_prefix + os.path.basename(infiles_paths[1]),
    ]
    for fn in outfiles:
        assert os.path.exists(fn)

    # Check that the number of duplicates is as expected
    stats_fn = "stats.csv"
    df_stats = pd.read_csv(stats_fn)
    n_duplicates = df_stats.query("stat_type == 'read_pairs_duplicated'").loc[:, "stat"].values[0]
    assert n_duplicates == n_duplicates_expected
