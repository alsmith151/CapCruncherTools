import os
from pathlib import Path

import pytest

from capcruncher_tools.api import deduplicate_fastq


@pytest.fixture(scope="module")
def data_path():
    return Path(__file__).resolve().parent / "fastq_deduplicate"


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
    data_path, tmp_path, monkeypatch, infiles, prefix, n_duplicates_expected
):
    monkeypatch.chdir(tmp_path)

    infiles_paths = [data_path / fn for fn in infiles]
    out_prefix = tmp_path / prefix

    df_stats = deduplicate_fastq(
        fastq1=[infiles_paths[0]],
        fastq2=[infiles_paths[1]],
        output_prefix=os.fspath(out_prefix),
    )

    # Check that the output files exist
    outfiles = [
        Path(f"{out_prefix}{os.path.basename(infiles_paths[0])}"),
        Path(f"{out_prefix}{os.path.basename(infiles_paths[1])}"),
    ]
    for fn in outfiles:
        assert fn.exists()

    # Check that the number of duplicates is as expected
    n_duplicates = df_stats.query("stat_type == 'read_pairs_duplicated'").loc[
        :, "stat"
    ].values[0]
    assert n_duplicates == n_duplicates_expected
