# Snakemake tutorial
This is the github repository for use with the session 
"Pipelining in Python: using Snakemake for biologicla applications". It shows 
how to use the workflow management system python. See the slides for more
information about how to use it.

## What is included
- Snakefile.smk includes the bulk of the code and is where the program is used.
This program uses utils.py, summarize_alignment.py, and config.yml.
- snakemake_tutorial_env.yml can be used to set up a conda enviornment with all
the needed programs. These programs can also be pip installed (instead of conda)
if that is your preference. More information below
- genomes and samples folders contain (relatively) small amounts of data needed 
to run the included code.

## Installing needed programs
I used conda, however it is possible to run everything without it 
(directly pip installing programs). If you use conda

#### create conda environment
```shell script
conda env create -f snakemake_tutorial_env.yml
```

#### update existing conda environment
```shell script
conda env update --file snakemake_tutorial_env.yml  --prune
```

#### update conda enviornment file: 
```shell script
conda env export --from-history > snakemake_tutorial_env.yml
```

## to cover
- submitting to cluster
- cluster.json
- saving dag

## Samples
- S5 bead bound
- S6 genomic dNA

### trouble shooting
if you get an error from samtools (mine was about libcrypto.1.0.0.dylib 
not being available), try running these three bash commands:
```shell script
conda uninstall samtools
conda update --all
conda install samtools
```