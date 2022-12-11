import logging
import os
import pathlib

import click
import pandas as pd

import capcruncher_tools.bindings

FORMAT = "%(levelname)s %(name)s %(asctime)-15s %(message)s"
logging.basicConfig(format=FORMAT)
logging.getLogger().setLevel(logging.DEBUG)


@click.group()
def cli():
    """CapCruncherTools CLI.
    Faster utilities to speed up CapCruncher
    """


@cli.command()
@click.option("-1", "--fastq1", help="Read 1 FASTQ files", required=True, multiple=True)
@click.option("-2", "--fastq2", help="Read 2 FASTQ files", required=True, multiple=True)
@click.option(
    "-o",
    "--output-prefix",
    help="Output prefix for deduplicated FASTQ files",
    default="deduped",
)
@click.option("-e", "--error-rate", help="Error rate for deduplication", default=0.01)
@click.option(
    "--shuffle", help="Shuffle shards before deduplication", default=False, is_flag=True
)
@click.option(
    "--sample-name", help="Name of sample e.g. DOX_treated_1", default="sampleX"
)
@click.option("--stats-prefix", help="Output prefix for stats file", default="stats")
def fastq_deduplicate(*args, **kwargs):

    fq_input = list(zip(kwargs["fastq1"], kwargs["fastq2"]))
    fq_output = [
        (
            os.path.join(kwargs["output_prefix"], os.path.basename(f1)),
            os.path.join(kwargs["output_prefix"], os.path.basename(f2)),
        )
        for f1, f2 in fq_input
    ]

    output_path = pathlib.Path(fq_output[0][0]).parent
    output_path.mkdir(parents=True, exist_ok=True)

    deduplication_results = capcruncher_tools.bindings.fastq_deduplicate(
        infiles=fq_input,
        outfiles=fq_output,
        shuffle=kwargs["shuffle"],
        error_rate=kwargs["error_rate"],
    )

    df_stats = (
        pd.Series(deduplication_results)
        .to_frame("stat")
        .reset_index()
        .rename(columns={"index": "stat_type"})
        .assign(
            read_number=0,
            read_type="pe",
            stage="deduplication",
            sample=kwargs["sample_name"],
        )
    )

    df_stats.to_csv(f"{kwargs['stats_prefix']}.deduplication.csv", index=False)

    print(df_stats)


if __name__ == "__main__":
    cli()
