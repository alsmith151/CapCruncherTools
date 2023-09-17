import os
import pathlib

import click
import pandas as pd
import tabulate
from loguru import logger as logging
from typing import Union, Tuple, List, Literal
import tempfile

from .capcruncher_tools import deduplicate, digest
from .count import count_interactions, make_cooler


import click
import ray


class OptionEatAll(click.Option):
    def __init__(self, *args, **kwargs):
        self.save_other_options = kwargs.pop("save_other_options", True)
        nargs = kwargs.pop("nargs", -1)
        assert nargs == -1, "nargs, if set, must be -1 not {}".format(nargs)
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
@click.option(
    "-1", "--fastq1", help="Read 1 FASTQ files", required=True, cls=OptionEatAll
)
@click.option(
    "-2", "--fastq2", help="Read 2 FASTQ files", required=True, cls=OptionEatAll
)
@click.option(
    "-o",
    "--output-prefix",
    help="Output prefix for deduplicated FASTQ files",
    default="deduped",
)
@click.option(
    "--sample-name", help="Name of sample e.g. DOX_treated_1", default="sampleX"
)
@click.option(
    "-s", "--statistics", help="Statistics output file name", default="stats.csv"
)
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

    logging.info(
        f"Digesting genome using recognition site: {kwargs['recognition_site']}"
    )
    logging.info(
        f"{'Keeping' if not kwargs['remove_recognition_site'] else 'Removing'} recognition site"
    )
    logging.info(f"Minimum slice length: {kwargs['min_slice_length']}")
    logging.info(f"Saving output to {kwargs['output']}")

    digest.digest_fasta(
        kwargs["input"],
        kwargs["recognition_site"],
        kwargs["output"],
        kwargs["remove_recognition_site"],
        kwargs["min_slice_length"],
        kwargs["n_threads"],
    )


@cli.command()
@click.argument("reporters")
@click.option("-o", "--output", help="Name of output file", default="CC_cooler.hdf5")
@click.option(
    "--remove_exclusions",
    default=False,
    help="Prevents analysis of fragments marked as proximity exclusions",
    is_flag=True,
)
@click.option(
    "--remove_capture",
    default=False,
    help="Prevents analysis of capture fragment interactions",
    is_flag=True,
)
@click.option(
    "--subsample",
    default=0,
    help="Subsamples reporters before analysis of interactions",
    type=float,
)
@click.option(
    "-f",
    "--fragment-map",
    help="Path to digested genome bed file",
)
@click.option(
    "-v",
    "--viewpoint-path",
    help="Path to viewpoints file",
)
@click.option(
    "-p",
    "--n-cores",
    default=1,
    help="Number of cores to use for counting.",
    type=int,
)
@click.option(
    "--assay", type=click.Choice(["capture", "tri", "tiled"]), default="capture"
)
def count(
    reporters: os.PathLike,
    output: os.PathLike = "counts.hdf5",
    remove_exclusions: bool = False,
    remove_viewpoint: bool = False,
    subsample: float = 0,
    fragment_map: os.PathLike = None,
    viewpoint_path: os.PathLike = None,
    n_cores: int = 1,
    assay: Literal["capture", "tri", "tiled"] = "capture",
    **kwargs,
):
    """Count interactions between restriction fragments in a supplied parquet file"""

    import pyranges as pr
    import capcruncher.api.storage
    import random
    import string
    import tqdm
    import tabulate

    df = pd.read_parquet(reporters, engine="pyarrow", columns=["viewpoint"])

    # Extract the viewpoint names and sizes from the parquet file
    viewpoints = df["viewpoint"].cat.categories.to_list()
    viewpoint_sizes = df["viewpoint"].value_counts()
    viewpoint_sizes_dict = viewpoint_sizes.to_dict()
    viewpoint_sizes_df = pd.DataFrame.from_dict(
        viewpoint_sizes_dict, orient="index", columns=["n_slices"]
    )
    viewpoint_sizes_df_tab = tabulate.tabulate(
        viewpoint_sizes_df, headers="keys", tablefmt="psql", showindex=True
    )


    logging.info(f"Number of viewpoints: {len(viewpoints)}")
    logging.info(f"Number of slices per viewpoint:\n {viewpoint_sizes_df_tab}")

    # If any viewpoint has more than 2 million slices, switch to low memory mode
    if any([vp for vp in viewpoint_sizes_dict.values() if vp > 2e6]):
        logging.warning(
            "High number of slices per viewpoint detected. Switching to low memory mode"
        )
        low_memory = True

        # Extract the partitions used to generate the file
        df = pd.read_parquet(reporters, engine="pyarrow", columns=["bam"])
        partitions = df["bam"].cat.categories.to_list()

    else:
        low_memory = False

    ray.init(num_cpus=n_cores, ignore_reinit_error=True)

    reporters_path = pathlib.Path(reporters)
    if reporters_path.is_dir():
        reporters = str(reporters_path / "*.parquet")

    counts = []
    for viewpoint in viewpoints:
        counts.append(
            count_interactions.remote(
                parquet=f"{reporters}",
                viewpoint=viewpoint,
                remove_exclusions=remove_exclusions,
                remove_viewpoint=remove_viewpoint,
                subsample=subsample,
                low_memory=low_memory,
                partitions=partitions if low_memory else None,
            )
        )

    bins = pr.read_bed(fragment_map, as_df=True).rename(
        columns={
            "Chromosome": "chrom",
            "Start": "start",
            "End": "end",
            "Name": "name",
        }
    )
    bins_ref = ray.put(bins)

    with tempfile.TemporaryDirectory() as tmpdir:
        cooler_refs = []
        for count_future in counts:
            cooler_refs.append(
                make_cooler.remote(
                    output_prefix=str(
                        pathlib.Path(tmpdir)
                        / f"{''.join(random.choices(string.ascii_letters + string.digits, k=8))}.hdf5"
                    ),
                    future=count_future,
                    bins=bins_ref,
                    viewpoint_path=viewpoint_path,
                    assay=assay,
                )
            )

        with tqdm.tqdm(total=len(cooler_refs)) as pbar:
            coolers = []
            while cooler_refs:
                coolers_completed, cooler_refs = ray.wait(cooler_refs)
                for clr in coolers_completed:
                    clr = ray.get(clr)
                    coolers.append(clr.split("::")[0])
                    pbar.update(1)

        logging.info(f"Making final cooler at {output}")
        capcruncher.api.storage.merge_coolers(coolers, output=output)


@cli.command("digest-fastq")
@click.option("-i","--input-fastqs", multiple=True, required=True)
@click.option(
    "-r",
    "--restriction_enzyme",
    help="Restriction enzyme name or sequence to use for in silico digestion.",
    required=True,
)
@click.option(
    "-m",
    "--mode",
    help="Digestion mode. Combined (Flashed) or non-combined (PE) read pairs.",
    type=click.Choice(["flashed", "pe"], case_sensitive=False),
    required=True,
)
@click.option("-o", "--output_file", default="out.fastq.gz")
@click.option("--minimum_slice_length", default=18, type=click.INT)
@click.option("--stats-prefix", help="Output prefix for stats file", default="stats")
@click.option(
    "--sample-name",
    help="Name of sample e.g. DOX_treated_1. Required for correct statistics.",
    default="sampleX",
)
def digest_fastq(*args, **kwargs):
    logging.info(f"Digesting FASTQ files using {kwargs['restriction_enzyme']}")
    (stats_read, stats_hist_unfilt, stats_hist_filt, stats_hist_len) = digest.digest_fastq(
        kwargs["input_fastqs"],
        kwargs["restriction_enzyme"].lower(),
        kwargs["output_file"],
        kwargs["mode"].capitalize(),
        kwargs["sample_name"],
        kwargs["minimum_slice_length"],
    )

    print("Read statistics")
    print(stats_read)
    stats_read.write_csv(kwargs["stats_prefix"] + ".digestion.csv")
    
    print("Histogram of unfiltered slice numbers")
    print(stats_hist_unfilt.sort(["slice_number", "count"]))
    
    print("Histogram of filtered slice numbers")
    print(stats_hist_filt.sort(["slice_number", "count"]))
    
    print("Histogram of filtered slice lengths")
    print(stats_hist_len.sort(["slice_length", "count"]))
    


if __name__ == "__main__":
    cli()
