# CapCruncherTools

## Overview

A collection of Rust tools to speed up the python based functionality of the CapCruncher project. Python binding have been generated to allow for easy integration into the existing python code base. 

## Current Tools

### FASTQ deduplication - fastq-deduplicate:

This tool takes paired FASTQ files and removes any duplicate fragments. 

#### Basic Usage:

```bash
capcruncher-tools fastq-deduplicate -1 <input1.fastq> -2 <input2.fastq> -o <output_prefix>
```


### Restriction digestion of FASTA - digest-genome:

This tool takes a FASTA file and a list of restriction enzymes and produces a list of fragments in BED format.

#### Basic Usage:

```bash
capcruncher-tools digest-genome -i <input.fasta> -o <output.bed> -r <recognition site> -p  <number of threads>
```


### Count restriction fragments - count:

This tool counts the number of interactions between a fragment (in silico digested read).

#### Basic Usage:

```bash
capcruncher-tools count <reporters> -f <fragments.bed> -v <viewpoints_path.bed> -o <output.hdf5> -p <number of threads>
```
