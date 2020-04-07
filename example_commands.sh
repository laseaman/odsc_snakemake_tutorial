# # setup conda environment for the first time
# conda env create -f snakemake_tutorial_env.yml

# # activate conda environment
conda activate snakemake_tutorial

# # dry run the pipeline
snakemake -s Snakefile.smk -n -r

# # run the pipeline
# this took just under 2 minutes to run on my 2019 macbook pro
snakemake -s Snakefile.smk --cores 1

# # saving an image of the rules to be run and their dependenceies (the Directed Acyclic Graph or DAG)
# this will fail if there are any print statements outside of rules. As a result, change verbose to
# False in the config before running it.
snakemake --forceall -s Snakefile.smk --dag | dot -Tpdf > dag.pd

# # run the pipeline - alternate version utalizing more of the command line options
# # --latency-wait 60: if at the end of any rule some of the files are missing, wait and check again
# #          after 60 seconds instead of immediately failing.
# # -k: keep running independant jobs after one fails (default is to submit no new jobs after a failure)
# # -R map_reads: fore rerun all map_read rule executions and all downstream jobs
# snakemake -s Snakefile.smk --cores 1 --latency-wait 60 -k -R map_reads

# # example of running jobs on a sbatch cluster
# snakemake --snakefile Snakefile.smk -j 12 --cluster "sbatch -n 1 -c 4 -p short --mem 4G -t 0-4:00

# # example of running jobs on a cluster using a json file to control the resources for each type of job
# snakemake --snakefile Snakefile.smk -j 12 --cluster-config cluster.json --cluster "sbatch -n {cluster.n} -c {cluster.c} -p short --mem {cluster.mem} -t {cluster.t}"