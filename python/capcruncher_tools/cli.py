import logging
import os
import pathlib

import click
import pandas as pd
import tabulate

from .capcruncher_tools import deduplicate, digest

FORMAT = "%(levelname)s %(name)s %(asctime)-15s %(message)s"
logging.basicConfig(format=FORMAT)
logging.getLogger().setLevel(logging.INFO)


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
    """Remove PCR duplicates from paired FASTQ files"""

    fq_input = list(zip(kwargs["fastq1"], kwargs["fastq2"]))
    fq_output = [
        (
            kwargs["output_prefix"] + os.path.basename(f1),
            kwargs["output_prefix"] + os.path.basename(f2),
        )
        for f1, f2 in fq_input
    ]

    output_path = pathlib.Path(fq_output[0][0]).parent
    output_path.mkdir(parents=True, exist_ok=True)

    deduplication_results = deduplicate.fastq_deduplicate(
        fq_input,  # Infiles
        fq_output,  # Outfiles
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

@cli.command()
@click.option("-i", "--input", help="Input fasta file", required=True)
@click.option("-o", "--output", help="Output bed file", default="digested.bed")
@click.option(
    "-r", "--recognition-site", help="Enzyme recognition site", default="GATC"
)
@click.option(
    "--min-slice-length",
    help="Minimum slice length to keep",
    default=18,
)
@click.option(
    "--remove-recognition-site",
    help="Remove recognition site from output",
    is_flag=True,
    default=False,
)
@click.option(
    "-p",
    "--n-threads",
    help="Number of threads to use for digesting",
    default=1,
)
def digest_genome(*args, **kwargs):
    """Digest genome with restriction enzyme and save to BED file"""

    # fasta: String,
    # restriction_site: String,
    # output: String,
    # remove_recognition_site: bool,
    # min_slice_length: Option<usize>,

    logging.info(f"Digesting genome using recognition site: {kwargs['recognition_site']}")
    logging.info(f"{'Keeping' if not kwargs['remove_recognition_site'] else 'Removing'} recognition site")
    logging.info(f"Minimum slice length: {kwargs['min_slice_length']}")
    logging.info(f"Saving output to {kwargs['output']}")
    
    digest.digest_fasta(
        kwargs["input"],
        kwargs["recognition_site"],
        kwargs["output"],
        kwargs["remove_recognition_site"],
        kwargs["min_slice_length"],
        kwargs["n_threads"]
    )

if __name__ == "__main__":
    cli()
