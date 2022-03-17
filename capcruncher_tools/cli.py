import logging

import click
import pandas as pd

import capcruncher_tools.bindings

FORMAT = "%(levelname)s %(name)s %(asctime)-15s %(message)s"
logging.basicConfig(format=FORMAT)
logging.getLogger().setLevel(logging.INFO)


@click.group()
def cli():
    """CapCruncherTools CLI.
    Faster utilities to speed up CapCruncher
    """


@cli.command()
@click.argument("infiles", nargs=2)
@click.option("-o", "--output-prefix", default="dd_")
@click.option("--output-compression", type=click.BOOL, default=False)
def fastq_deduplicate(*args, **kwargs):

    deduplication_results = capcruncher_tools.bindings.fastq_deduplicate(
        *args, **kwargs
    )

    logging.info(pd.Series(deduplication_results).to_frame("count"))


if __name__ == "__main__":
    cli()
