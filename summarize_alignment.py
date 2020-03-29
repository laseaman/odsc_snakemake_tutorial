#!/usr/bin/env python3
###########################################################################
# Purpose of this program is to tabulate alignment statistics for multiple files.
#   it makes a table of basic alignment statistics like overall alignment rate, number
#   of reads, and number of half mapped reads. It also makes tables of the total
#   number and RMP for each chromosome/contig/probe.
###########################################################################

# Libraries
import os
from pathlib import Path
import argparse
from warnings import warn
from io import StringIO
import pandas as pd
import sys
import re

BOWTIE2_COLS = ['num_reads', 'percent_concordantly_0times', 'percent_concordantly_1time',
                'percent_concordantly_2plus_times', 'pairs_concordantly_0times', 'pairs_discordantly_1time',
                'pairs_unaligned', 'reads_unaligned', 'reads_halfmapped', 'overall_alignment_rate']
SAMTOOLS_IDXSTATS_COLS = ['chr', 'num_reads', 'ignore']


def find_sample_id(file_name: str,
                   pattern: str) -> str:
    """Gets the sample id for the given file.
    :param file_name: File name to get sample id from.
    :param pattern: Glob pattern which matches the file_name, where wildcards match to the sample id.
      It is currently assumed that there is a single instance of the sample id in the file_name, and
      therefore a single wildcard in pattern.
    :example: find_sample_id('3H.align.o', '*.align.o') --> '3H'
    """
    n_id_instances = pattern.count('*')
    if n_id_instances == 0:
        raise ValueError('Glob pattern must have a wildcard to specify location of sample ID.')
    elif n_id_instances > 1:
        warn('Multiple wildcards in glob pattern; may cause unexpected behavior')

    sample_id = file_name
    for s in pattern.split('*'):
        sample_id = sample_id.replace(s, '', 1)
    return sample_id


def parse_bowtie2_summary(text, columns):
    """
    This will parse the .o (std out) file produced during alignment and return the interesting
        information in a pandas series.
    :param text: the complete text of the alignment output file
    :param columns: list of columns to be included in the dataframe
    :return: pandas dataframe of results with a single row for a single sample
    """
    # TODO: print a message if a column is empty because it does not have an if statement below
    results = pd.DataFrame(columns=columns)
    if len(text) < 10:
        print('It seems like alignment failed or had not completed becuase the alignment summary is not complete.')
        print(text)
        data = pd.DataFrame(columns=columns, data=0, index=[0])
        return data
    # python command used for alignment
    if 'num_reads' in columns:
        line = text[text.str.contains('reads; of these:')]
        results.loc[0, 'num_reads'] = int(line.iloc[0].split()[0])
    # percent concordantly 0 times (unmapped or discordant, i.e. wrong distance or strand)
    if 'percent_concordantly_0times' in columns:
        line = text[text.str.contains('aligned concordantly 0 times')]
        val = line.iloc[0].split(')')[0].split('(')[-1]
        results.loc[0, 'percent_concordantly_0times'] = float(val.strip('%'))
    if 'percent_concordantly_1time' in columns:
        line = text[text.str.contains('aligned concordantly exactly 1 time')]
        val = line.iloc[0].split(')')[0].split('(')[-1]
        results.loc[0, 'percent_concordantly_1time'] = float(val.strip('%'))
    if 'percent_concordantly_2plus_times' in columns:
        line = text[text.str.contains('aligned concordantly >1 times')]
        val = line.iloc[0].split(')')[0].split('(')[-1]
        results.loc[0, 'percent_concordantly_2plus_times'] = float(val.strip('%'))
    if 'pairs_concordantly_0times' in columns:
        line = text[text.str.contains('aligned concordantly 0 times; of these:')]
        results.loc[0, 'pairs_concordantly_0times'] = int(line.iloc[0].split()[0])
    if 'pairs_discordantly_1time' in columns:
        line = text[text.str.contains('aligned discordantly 1 time')]
        results.loc[0, 'pairs_discordantly_1time'] = int(line.iloc[0].split()[0])
    # read pairs that did not align (at least one read)
    if 'pairs_unaligned' in columns:
        line = text[text.str.contains('pairs aligned 0 times concordantly or discordantly')]
        results.loc[0, 'pairs_unaligned'] = int(line.iloc[0].split()[0])
    if 'reads_unaligned' in columns:
        line = text[text.str.contains('aligned 0 times')]
        results.loc[0, 'reads_unaligned'] = int(line.iloc[1].split()[0])
    # the number of individual reads that did align from read pairs that did not align
    # this is the number of half-mapped reads, because their mate did not map, but they did
    if 'reads_halfmapped' in columns:
        line = text[text.str.contains('aligned exactly 1 time')]
        line2 = text[text.str.contains('aligned >1 times')]
        results.loc[0, 'reads_halfmapped'] = int(line.iloc[0].split()[0]) + int(line2.iloc[0].split()[0])
    if 'overall_alignment_rate' in columns:
        line = text[text.str.contains('overall alignment rate')]
        val = line.iloc[0].split()[0]
        results.loc[0, 'overall_alignment_rate'] = float(val.strip('%'))
    return results


def run_return(cmd, to_table=True, columns=None):
    """
    This will run a command in a way that returns the output to python so it can be analyzed.
    :param cmd: list of string to be joined with spaces that will form the final command
    :param to_table: binary for if the results are a table and therefor should be converted into a pandas dataframe
    :param columns: the column names for the resulting data frame.
    :return: the contents of the command as a string (if to_table=False) or a pandas dataframe (if to_table=True)
    """
    string = os.popen(' '.join(cmd)).read()
    if to_table:
        table = pd.read_csv(StringIO(string), names=columns, sep='\t')
        return table
    return string.strip()


if __name__ == '__main__':
    # define parser and parse args
    parser = argparse.ArgumentParser(description='Tabulates alignment statistics for multiple files.')
    # Required Arguments
    parser.add_argument('alignment_dir',
                        help='path to target input directory containing alignment files',)
    parser.add_argument('log_dir',
                        help='path to target input directory containing alignment log files', )
    parser.add_argument('out_dir',
                        help='output directory')
    # Optional Arguments
    parser.add_argument('-d', '--description',
                        help='description of the alignment, or the name of the run; used to name output files')
    parser.add_argument('-lf', '--logfile_pattern', default='*align_to_probes_stats.log',
                        help='glob pattern matching the alignment statistics files')
    parser.add_argument('-l', '--logfile', default=None,
                        help='the file name to which print statements should be directed. Use none to print to stdout')
    parser.add_argument('-b', '--bam_pattern', default='_sorted.bam',
                        help='glob pattern matching the sorted .bam files')
    parser.add_argument('-v', '--verbose', default='False',
                        help='verbose = True will print out things to help with debugging')
    parser.add_argument('-p', '--probe_ref_genome', default=None,
                        help='will get length of probes used in run for RPKM table')
    args = parser.parse_args()

    # Initialization
    alignment_dir = Path(args.alignment_dir)
    log_dir = Path(args.log_dir)
    out_dir = Path(args.out_dir)
    alignstats_outfile = 'alignment_stats.csv'
    reads_outfile = 'reads_per_chr.csv'
    rpm_outfile = 'rpm_per_chr.csv'
    rpkm_outfile = 'rpkm_per_chr.csv'
    if args.description is not None:
        alignstats_outfile = '_'.join([args.description, alignstats_outfile])
        reads_outfile = '_'.join([args.description, reads_outfile])
        rpm_outfile = '_'.join([args.description, rpm_outfile])
        rpkm_outfile = '_'.join([args.description, rpkm_outfile])
    # open the log file
    if args.logfile is not None:
        log_file = open(args.logfile, "w")
        sys.stdout = log_file

    ''' PART 1: read through alignment files to get statistics '''
    align_data = pd.DataFrame(columns=BOWTIE2_COLS)
    if args.verbose:
        print(args.logfile_pattern)
    for f in list(log_dir.glob(args.logfile_pattern)):
        sample_id = find_sample_id(f.name, args.logfile_pattern)
        with f.open('r') as fp:
            lines = pd.Series(fp.readlines())
            tmp = parse_bowtie2_summary(lines, BOWTIE2_COLS)
            align_data.loc[sample_id] = tmp.loc[0]
    align_data.sort_index(inplace=True)
    out_path = out_dir / alignstats_outfile
    if args.verbose:
        print(align_data.iloc[0:3, 2:])
        print(f'saving to {out_path}')
    align_data.to_csv(out_path)

    ''' PART 2: read through bam files to find distribution of reads '''
    # set up the dataframe for collecting results
    chr_table = pd.DataFrame(columns=['sample', 'total_reads'])
    # count number of reads aligning to each probe
    bam_files = list(alignment_dir.glob(args.bam_pattern))
    if args.verbose:
        print(bam_files)
        print(alignment_dir)
        print(args.bam_pattern)
    for i in range(len(bam_files)):
        f = bam_files[i]
        sample_id = find_sample_id(f.name, args.bam_pattern)
        # count reads per chr
        reads_table = run_return(['samtools idxstats ', str(f)], columns=SAMTOOLS_IDXSTATS_COLS)
        # convert format and put in table
        chr_table = chr_table.append(reads_table.num_reads, ignore_index=True)
        chr_table.loc[i, 'sample'] = sample_id
        # get total number of reads
        num_reads = run_return(['samtools view -c ', str(f)], to_table=False)
        chr_table.loc[i, 'total_reads'] = int(num_reads)

    # save the raw count table
    table = chr_table.set_index('sample').sort_index()
    out_path = out_dir / reads_outfile
    if args.verbose:
        print(table.head())
        print(f'saving to {out_path}')
    table.to_csv(out_path)

    # create a table of normalized counts (reads per million = RPM)
    norm_table = table.div(table.total_reads, axis=0) * 1000000
    norm_table.total_reads = table.total_reads
    out_path = out_dir / rpm_outfile
    if args.verbose:
        print(f'saving to {out_path}')
    norm_table.to_csv(out_path)

    # create a table for dividing by length of each probe and multiply by 1000
    if args.probe_ref_genome is not None:
        probe_file_path = args.probe_ref_genome + '.fa'
        os.system('samtools faidx ' + probe_file_path)
        os.system('cut -f1,2 ' + probe_file_path + '.fai > results/sizes_genome_' + args.description + '.txt')
        with open("results/sizes_genome_" + args.description + ".txt", 'r') as f:
            size_list = f.read().splitlines()
        size_list.sort()
        numbers_only = [1, 1]
        for element in size_list:
            numbers_only.append(int(re.findall("[0-9]*$", element)[0]))
        if args.verbose:
            print('Dividing columns by: ')
            print(numbers_only)
        RPKM_table = table.div(numbers_only, axis=1) * 1000
        RPKM_table.total_reads = table.total_reads
        out_path = out_dir / rpkm_outfile
        if args.verbose:
            print(f'saving to {out_path}')
        RPKM_table.to_csv(out_path)

    # close the log file
    if args.logfile is not None:
        sys.stdout = sys.__stdout__
        log_file.close()
