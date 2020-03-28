# ODSC 2020 snakemake tutorial

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

pip install needed programs
TODO: Add directions

update conda enviornment file: 
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