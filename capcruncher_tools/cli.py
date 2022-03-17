import logging
import pathlib

import click
import pandas as pd

from . import cct
import capcruncher_tools.bindings

FORMAT = '%(levelname)s %(name)s %(asctime)-15s %(filename)s:%(lineno)d %(message)s'
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
def fastq_deduplicate(*args, **kwargs):
    
    deduplication_results = capcruncher_tools.bindings.fastq_deduplicate(
        *args, **kwargs
    )

    print(pd.Series(deduplication_results))


if __name__ == "__main__":
    cli()
