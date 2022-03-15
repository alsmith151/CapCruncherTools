import click

@click.group()
def cli():
    """CapCruncherTools CLI. 
    Faster utilities to speed up CapCruncher
    """


@cli.command()
@click.argument("infiles", nargs=2)
@click.option("-o", "--output-prefix", default="dd_")
def fastq_deduplicate(infiles: tuple, output_prefix: str):
    """Deduplicates paired FASTQ files.
       Input files in format: fq1_1,fq2_1,fq3_1 fq1_2,fq2_2,fq3_2 
    """

    import capcruncher_tools

    fq_in = {
        str(ii): [r1, r2]
        for ii, (r1, r2) in enumerate(zip(*[files.split(",") for files in infiles]))
    }

    fq_out = {
        str(ii): [f"{output_prefix}{read}" for read in (r1, r2)]
        for ii, (r1, r2) in enumerate(zip(*[files.split(",") for files in infiles]))
    }

    capcruncher_tools.deduplicate_fastq(fq_in, fq_out)



if __name__ == "__main__":
    cli()