# Snakemake tutorial
This is the github repository for use with the session 
"Pipelining in Python: using Snakemake for biologicla applications". It shows 
how to use the workflow management system python. See the slides for more
information about how to use it.

More documentation and tutorials are available on the snakemake website: 
https://snakemake.readthedocs.io/en/stable/ 

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

create conda environment
```shell script
conda env create -f snakemake_tutorial_env.yml
```

update existing conda environment
```shell script
conda env update --file snakemake_tutorial_env.yml  --prune
```

update conda enviornment file: 
```shell script
conda env export --from-history > snakemake_tutorial_env.yml
```

## runing the code
- example commands including a dry run, full run, use of command line parameters, 
and cluster submission is in example_commands.sh
- to run the snakemake pipeline, first make sure the dependencies are installed and/or
you have activated the conda environment then run:
```shell script
snakemake -s Snakefile.smk --cores 1
```

## Included data
- pseudomonas.fa: a text file (fasta format), containg the genetic sequnce for 
pseudomonas aeruginosa a bacteria that is commonly used for genetic studies.
- S5_bead_L001_R1_001.fsatq.gz - gzipped fastq (text file for storing sequencind data) file
- S5_bead_L001_R1_001.fsatq.gz - same as S5_bead_L001_R1_001.fsatq.gz but the data for read 2 in the
     paired reads (background information: https://systemsbiology.columbia.edu/genome-sequencing-defining-your-experiment)
- S6_control_L001_R1_001.fsatq.gz - same as above but for a control sample (read 1)
- S6_control_L001_R1_001.fsatq.gz - same as above but for a control sample (read 2)

### trouble shooting
if you get an error from samtools (mine was about libcrypto.1.0.0.dylib 
not being available), try running these three bash commands:
```shell script
conda uninstall samtools
conda update --all
conda install samtools
```
