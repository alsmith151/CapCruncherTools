import logging
import os
import pathlib

import click
import pandas as pd
import tabulate

from .capcruncher_tools import deduplicate

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
@click.option(
    "--sample-name", help="Name of sample e.g. DOX_treated_1", default="sampleX"
)
@click.option("--stats-prefix", help="Output prefix for stats file", default="stats")
@click.option(
    "--shuffle",
    help="Shuffle reads before deduplication",
    is_flag=True,
    default=False,
)
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

    deduplication_results = deduplicate.fastq_deduplicate(
        fq_input, # Infiles
        fq_output, # Outfiles
        kwargs["shuffle"],
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

    logging.info(f"Saving stats to {kwargs['stats_prefix']}.deduplication.csv")
    df_stats.to_csv(f"{kwargs['stats_prefix']}.deduplication.csv", index=False)

    logging.info("Printing deduplication statistics to stdout")
    # Print stats to stdout
    df_vis = df_stats.copy()
    df_vis["stat_type"] = df_vis["stat_type"].str.replace("_", " ").str.title()
    df_vis = df_vis[["stat_type", "stat"]]
    df_vis.columns = ["Stat Type", "Number of Reads"]
    print(tabulate.tabulate(df_vis, headers="keys", tablefmt="psql", showindex=False))


if __name__ == "__main__":
    cli()
