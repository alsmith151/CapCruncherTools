# CapCruncherTools

## Overview

A collection of Rust tools to speed up the python based functionality of the CapCruncher project. Python binding have been generated to allow for easy integration into the existing python code base. 

## Current Accelerators

### FASTQ deduplication

This tool takes paired FASTQ files and removes any duplicate fragments. 

Use it through the CapCruncher CLI or Python API; this package only provides the
compiled acceleration layer.


### Restriction digestion of FASTA

This tool takes a FASTA file and a list of restriction enzymes and produces a list of fragments in BED format.

Use it through the CapCruncher CLI or Python API; this package does not install a
separate `capcruncher-tools` command.
