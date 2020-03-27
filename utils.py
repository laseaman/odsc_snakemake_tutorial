import re
from itertools import chain


def load_samples(config, verbose):
    # initiate things
    samples_dict = config["SAMPLES"]
    # the key is what to use when saving results
    description = list(samples_dict.keys())
    orig_name = list(samples_dict.values())
    # create list of paths to raw data
    all_fastq =list(chain.from_iterable((x + '1_001', x + '2_001') for x in orig_name))

    if verbose:
        print('dictionary from config: ', samples_dict)
        print('descriptions to use when saving: ', description)
        print('original file names: ', all_fastq)
    return description, all_fastq