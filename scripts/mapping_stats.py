#!/usr/bin/env python3
#usage: python3 scripts/stats_mapping.py -i sample_test.tsv -o test_out.tsv

# imports:
import os
import argparse
import pandas as pd
import numpy as np
import re
import sys

def read_file(file):
    """ Read the sample file containing path files and return a pandas dataframe"""
    sample = pd.read_csv(file, sep='\t', dtype=str)
    return sample

def get_log_stats(log_file):
    """ Reading, parsing the log file and returning a dict"""
    log_statistics = []
    with open(log_file) as fin:
        for line in fin:
            line = line.lstrip().rstrip()
            if not line or line.endswith(":") :  # if the line is empty we simply skip it
                continue
            key, value = line.split("|")
            value = value.replace("%", "")
            stat = (key.lstrip().rstrip(), value.lstrip().rstrip())
            log_statistics.append(stat)
    return log_statistics


def get_mapping_stats(runs):
    """
        Parse the STAR Log.final.out for each run to gather alignments
        statistics and returns a pandas dataframe
    """
    statistics = []
    reads = ["R1", "R2"]
    for idx, row in runs.iterrows():
        sample_info = [('sample', row['sample']),
                        ('sample_unit',row['sample_unit'] )]
        for read in reads:
            mapdir = os.path.join(row["mapdir"], read)
            if not os.path.exists(mapdir):
                raise Exception("WARNING mapdir %s is missing" % mapdir)
            log_file = os.path.join(mapdir, "Log.final.out")
            if not os.path.exists(log_file):
                raise Exception("WARNING logfile %s is missing" % log_file)
            sample_info.append(('read', read ))
            log_stats = get_log_stats(log_file)
            sample_stats = sample_info + log_stats

            statistics.append(dict(sample_stats))
    return pd.DataFrame(statistics)


def main(input, output):

    # Read the sample file:
    runs = read_file(input)

    # Reading the log files
    stats = get_mapping_stats(runs)

    # Writing statistics to a file in tsv format
    stats.to_csv(output, sep="\t", index=False)



def parse_arguments():
    parser = argparse.ArgumentParser(description='Sample file')
    parser.add_argument('-i', '--input_file',
                        required=True, help='Sample file')
    parser.add_argument('-o', '--output_file',
                        required=True, help='Statuistics in tsv format')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main(args.input_file, args.output_file)
