import re
from itertools import chain


def load_samples(config, verbose):
    # initiate things
    samples = config["SAMPLES"]
    # create list of paths to raw data
    r1 = config["fastq_suffix"][0]
    r2 = config["fastq_suffix"][1]
    all_fastq = list(chain.from_iterable((x + r1, x + r2) for x in samples))

    if verbose:
        print('the config file as loaded: ')
        print(config)
        print('descriptions to use when saving: ', samples)
        print('full file names: ', all_fastq)
    return samples, all_fastq