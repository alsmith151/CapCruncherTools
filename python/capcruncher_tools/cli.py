import os
import pathlib

import click
import pandas as pd
import tabulate
from loguru import logger as logging

from .capcruncher_tools import deduplicate, digest

import click

class OptionEatAll(click.Option):

    def __init__(self, *args, **kwargs):
        self.save_other_options = kwargs.pop('save_other_options', True)
        nargs = kwargs.pop('nargs', -1)
        assert nargs == -1, 'nargs, if set, must be -1 not {}'.format(nargs)
        super(OptionEatAll, self).__init__(*args, **kwargs)
        self._previous_parser_process = None
        self._eat_all_parser = None

    def add_to_parser(self, parser, ctx):

        def parser_process(value, state):
            # method to hook to the parser.process
            done = False
            value = [value]
            if self.save_other_options:
                # grab everything up to the next option
                while state.rargs and not done:
                    for prefix in self._eat_all_parser.prefixes:
                        if state.rargs[0].startswith(prefix):
                            done = True
                    if not done:
                        value.append(state.rargs.pop(0))
            else:
                # grab everything remaining
                value += state.rargs
                state.rargs[:] = []
            value = tuple(value)

            # call the actual process
            self._previous_parser_process(value, state)

        retval = super(OptionEatAll, self).add_to_parser(parser, ctx)
        for name in self.opts:
            our_parser = parser._long_opt.get(name) or parser._short_opt.get(name)
            if our_parser:
                self._eat_all_parser = our_parser
                self._previous_parser_process = our_parser.process
                our_parser.process = parser_process
                break
        return retval


@click.group()
def cli():
    """CapCruncherTools CLI.
    Faster utilities to speed up CapCruncher
    """


@cli.command()
@click.option("-1", "--fastq1", help="Read 1 FASTQ files", required=True, cls=OptionEatAll)
@click.option("-2", "--fastq2", help="Read 2 FASTQ files", required=True, cls=OptionEatAll)
@click.option(
    "-o",
    "--output-prefix",
    help="Output prefix for deduplicated FASTQ files",
    default="deduped",
)
@click.option(
    "--sample-name", help="Name of sample e.g. DOX_treated_1", default="sampleX"
)
@click.option("-s", "--statistics", help="Statistics output file name", default="stats.csv")
@click.option(
    "--shuffle",
    help="Shuffle reads before deduplication",
    is_flag=True,
    default=False,
)
def fastq_deduplicate(*args, **kwargs):
    """Remove PCR duplicates from paired FASTQ files"""

    logging.info(f"Output prefix: {kwargs['output_prefix']}")
    logging.info(f"Sample name: {kwargs['sample_name']}")
    logging.info(f"Stats prefix: {kwargs['statistics']}")
    logging.info(f"Shuffle reads: {kwargs['shuffle']}")

    import ast
    
    fq1 = ast.literal_eval(kwargs["fastq1"])
    fq2 = ast.literal_eval(kwargs["fastq2"])
        
    fq_input = list(zip(fq1, fq2))

    fq_output = [
        (
            kwargs["output_prefix"] + os.path.basename(f1),
            kwargs["output_prefix"] + os.path.basename(f2),
        )
        for f1, f2 in fq_input
    ]

    output_path = pathlib.Path(fq_output[0][0]).parent
    output_path.mkdir(parents=True, exist_ok=True)

    stats_path = pathlib.Path(kwargs["statistics"])

    if not stats_path.parent.exists():
        raise ValueError(f"Statistics path {stats_path.parent} does not exist")
    

    logging.info("Running deduplication")
    deduplication_results = deduplicate.fastq_deduplicate(
        fq_input,  # Infiles
        fq_output,  # Outfiles
        kwargs["shuffle"],
    )

    logging.info("Saving deduplication statistics")
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

    logging.info(f"Saving stats to {stats_path}.deduplication.csv")
    df_stats.to_csv(stats_path, index=False)

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
