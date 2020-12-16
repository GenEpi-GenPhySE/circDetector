#!/usr/bin/env python3
#usage: python3 scripts/stats_mapping.py -i samples.tsv -o mapping_stats_samples.tsv

# imports:
import os
from os import path
import argparse
import pandas as pd
import numpy as np
import re
import sys

from collections import defaultdict


def eprint(*args, **kwargs):
    print(*args,  file=sys.stderr, **kwargs)


def main(bsj_counts, fsj_counts, output_file):

    bsj = pd.read_csv(bsj_counts, sep='\t').set_index('bsj').drop(['gene_id', 'circ_type'], axis=1)
    fsj = pd.read_csv(fsj_counts, sep='\t').set_index('fsj').drop(['gene_id', 'circ_type'], axis=1)

    merged = pd.merge(bsj, fsj, how='inner', left_index=True, right_index=True,
                      suffixes=['_bsj', '_fsj'])

    merged.to_csv(output_file, sep="\t")


def parse_arguments():
    parser = argparse.ArgumentParser(description='Merge bsj and fsj counts')
    parser.add_argument('-b', '--bsj',
                        required=True, help='bsj count file')
    parser.add_argument('-f', '--fsj',
                        required=True, help='fsj count file')
    parser.add_argument('-o', '--output_file',
                        required=True, help='output file')


    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main(args.bsj, args.fsj, args.output_file)
