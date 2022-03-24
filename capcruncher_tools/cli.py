from builtins import breakpoint
import logging

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
@click.argument("infiles", nargs=2)
@click.option("-o", "--output-prefix", default="dd_")
@click.option("--sample-name", help="Name of sample e.g. DOX_treated_1", default='sampleX')
@click.option("--stats-prefix", help="Output prefix for stats file", default='stats')
@click.option("--output-compression", type=click.BOOL, default=False)
def fastq_deduplicate(*args, **kwargs):

    deduplication_results = capcruncher_tools.bindings.fastq_deduplicate(
        *args, **kwargs
    )

    df_stats = (pd.Series(deduplication_results)
                 .to_frame("stat")
                 .reset_index()
                 .rename(columns={"index": "stat_type"})
                 .assign(read_number=0,
                         read_type="pe",
                         stage="deduplication",
                         sample=kwargs["sample_name"], 
                         ))
    
    df_stats.to_csv(f"{kwargs['stats_prefix']}.deduplication.csv", index=False)


if __name__ == "__main__":
    cli()
