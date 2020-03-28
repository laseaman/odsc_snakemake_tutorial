# Purpose:  performs sequence alignment
# Run:      conda activate snakemake-tutorial; snakemake -s Snakefile.smk
# Dry-Run:      conda activate snakemake-env; snakemake -s Snakefile.smk -n
# Submission Example:
# snakemake --snakefile Snakefile.smk -j 12 --cluster "sbatch -n 4 -c 1 -p short --mem 4G -t 0-4:00" --latency-wait 60

import utils

##-----------------------------------------------##
## Config and Helper Files                       ##
##-----------------------------------------------##
configfile: "config.yml"
print(config)
##-----------------------------------------------##
## Access Variables from Config File             ##
##-----------------------------------------------##
# unpack the config file for later use
BASE_DIR = config["BASE_DIR"]
workdir: BASE_DIR
REF_DIR = config["REF_DIR"]
GENOME_NAME = config["REF"]
SAVE_LOC = config['SAVE_LOC']
core = config["CORES"]

# Get sample names (e.g. "S5")
verbose = config['verbose']
IDS, full_names = utils.load_samples(config, verbose=verbose)

# TODO: add visualization if time


rule all:
    input:
        # check quality
       expand("fastqc/{full_name}_fastqc.zip", full_name = full_names),
        # alignment stats summary
        "results/probes_alignment_stats.csv",

rule run_fastqc:
    input:
       expand("samples/{full_name}"
              ".fastq.gz", full_name = full_names)
    output:
        expand("fastqc/{full_name}_fastqc.html", full_name = full_names),
        expand("fastqc/{full_name}_fastqc.zip", full_name = full_names)
    log:
        expand("log_files/fastqc/{full_name}.log", full_name = full_names)
    params:
        folder=directory("fastqc/")
    shell:
        """
        for file in {input}
        do
            fastqc --outdir={params.folder} $file
        done
        """

rule make_bowtie2_index:
    output:
        bt2_index = '{REF_DIR}{REF_GENOME}.1.bt2'
    log:
        'log_files/{REF_DIR}{REF_GENOME}_make_bt2.log'
    shell:
        """
        (bowtie2-build {REF_DIR}/{GENOME_NAME}.fa {REF_DIR}/{GENOME_NAME}) &> {log}
        """


rule map_reads:
    input:
        R1="samples/{sample}_L001_R1_001.fastq.gz",
        R2="samples/{sample}_L001_R2_001.fastq.gz",
        bt2_index = expand('{ref_dir}/{ref}.1.bt2', ref_dir=REF_DIR, ref=GENOME_NAME)
    output:
        sam=temp("aligned_reads/{sample}_aligned.sam"),
        bam = "aligned_reads/{sample}_aligned.bam",
        sorted = 'aligned_reads/{sample}_aligned_sorted.bam',
    log:
        "log_files/{sample}_aligned_stats.log"
    shell:
        """
        (bowtie2 -x {REF_DIR}/{GENOME_NAME} -p {core} -1 {input.R1} -2 {input.R2} -S {output.sam}) 2> {log}
        samtools view -bS {output.sam} > {output.bam}
        samtools sort {output.bam} -o {output.sorted}
        """


# parse alignment stats to create a table with the percent aligned and a table with the number of reads mapping to each
#    engineering maker
rule parse_alignment_stats:
    input:
        bams=expand('aligned_reads/{sample}_aligned_sorted.bam', sample=IDS),
    params:
        bam_dir=directory("aligned_reads/"),
        log_dir=directory("log_files/"),
        extension='*aligned_stats.log',
        extension2='*_aligned_sorted.bam',
        out_fold='results/',
    log:
        "log_files/parse_alignment_stats.log"
    output:
        stat="results/probes_alignment_stats.csv",
        read="results/probes_reads_per_chr.csv",
        rpm="results/probes_rpm_per_chr.csv",
        rpkm="results/probes_rpkm_per_chr.csv"
    shell:
        """
        (python summarize_alignment.py {params.bam_dir} {params.log_dir} {params.out_fold} -d probes -lf {params.extension} -b {params.extension2} -v {verbose} -p {REF_DIR}/{GENOME_NAME}) &>{log}
        """
