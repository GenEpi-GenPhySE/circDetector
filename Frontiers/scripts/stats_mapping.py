#!/usr/bin/env python3
#usage: python3 scripts/stats_mapping.py -i samples.tsv -o mapping_stats_samples.tsv

# imports:
import os
import argparse
import pandas as pd
import numpy as np
import re
import sys

from collections import defaultdict


# Utility functions
def create_folder(path):
    if not os.path.exists(path):
        os.mkdir(path)

def eprint(*args, **kwargs):
    print(*args,  file=sys.stderr, **kwargs)


def set_stats_dictionnary():
    d = {
    	'Started job on ': 'start',
        'Started mapping on ': 'start mapping',
        'Finished on ': 'end',
        'Mapping speed, Million of reads per hour ': 'speed',
        'Number of input reads ': 'nb input reads',
        'Average input read length ': 'readlen',
        'Uniquely mapped reads number ': 'uniquely mapped',
        'Uniquely mapped reads % ': 'uniquely mapped %',
        'Average mapped length ': 'average map length',
        'Number of splices: Total ': 'nb splices',
        'Number of splices: Annotated (sjdb) ': 'nb annotated splices',
        'Number of splices: GT/AG ': 'nb splices GT/AG',
        'Number of splices: GC/AG ': 'nb splices GC/AG',
        'Number of splices: AT/AC ': 'nb splices AT/AC',
        'Number of splices: Non-canonical ': 'nb splices noncanon',
        'Mismatch rate per base, % ': 'mismatch rate',
        'Deletion rate per base ': 'deletion rate',
        'Deletion average length ': 'del avg len',
        'Insertion rate per base ': 'ins rate',
        'Insertion average length ': 'ins avg len',
        'Number of reads mapped to multiple loci ': 'multiple loci',
        '% of reads mapped to multiple loci ': 'multiple loci %',
        'Number of reads mapped to too many loci ': 'too many loci',
        '% of reads mapped to too many loci ': 'too many loci %',
        '% of reads unmapped: too many mismatches ': 'unmapped, too many mismatches % ',
        '% of reads unmapped: too short ': 'unmapped, too short %',
        '% of reads unmapped: other ': 'unmapped, other %',
        'Number of chimeric reads ': 'chimeric reads',
        '% of chimeric reads ': 'chimeric reads %'
	}
    return d

def get_valid_keys():
    keys = ['nb input reads', 'readlen', 'uniquely mapped', 'uniquely mapped %',
            'multiple loci', 'multiple loci %', 'too many loci',
            'too many loci %', 'chimeric reads']
    return keys


def read_file(file):
    """ Read the sample file containing path files and return a pandas dataframe"""
    sample = pd.read_csv(file, sep='\t', dtype=str)
    return sample


def sample_stats(sample, species, root_mapping_dir):
    mapping_dir = os.path.join(root_mapping_dir, species, "mapping", sample)

    stats = defaultdict()
    # se
    for orient in ["R1", "R2"]:
        logfile = os.path.join(mapping_dir, "se", orient, "Log.final.out")
        stats[orient] = parse_log_file(logfile)
    logfile = os.path.join(mapping_dir, "pe", "Log.final.out")
    stats["pe"] = parse_log_files(logfile)
    return stats


def parse_log_file(logfile):
    stats = defaultdict()
    d_keys = set_stats_dictionnary()
    valid_keys = get_valid_keys()
    with open(logfile, "r") as f:
        for line in f:
            fields = line.rstrip().lstrip().split("|\t")
            if len(fields) == 1:
                continue
            if fields[0] in d_keys and d_keys[fields[0]] in valid_keys:
                stats[d_keys[fields[0]]] = fields[1]
    return stats


def stats_to_dataframe(stats):
    valid_keys = get_valid_keys()
    array = []
    for sample in sorted(stats):
        for type, values in stats[sample].items():
            infos = [sample, type] + [values[k] for k in valid_keys]
            array.append(infos)
    columns = ["sample", "type"] + valid_keys
    df = pd.DataFrame(array, columns=columns)
    return df

def main(input_file, input_dir, output_file):

    samples = pd.read_csv(input_file, sep='\t', dtype=str).set_index("sample")

    stats = defaultdict()
    for index, row in samples.iterrows():
        stats[index] = sample_stats(index, row["species"], input_dir)

    df = stats_to_dataframe(stats)
    df.to_csv(output_file, sep="\t", index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Sample file')
    parser.add_argument('-i', '--input_file',
                        required=True, help='Sample file')
    parser.add_argument('-d', '--input_dir',
                        required=True, help='root dir of the mapping results')
    parser.add_argument('-o', '--output_file', required=False,
                        default="mapping_stat.tsv", help='Sample file')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main(args.input_file, args.input_dir, args.output_file)
