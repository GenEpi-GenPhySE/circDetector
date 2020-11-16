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

def eprint(*args, **kwargs):
    print(*args,  file=sys.stderr, **kwargs)

def main(input_file, output_file):

    circrnas = pd.read_csv(input_file, sep='\t')

    counts4 = (circrnas.drop(['circ','genename'], axis=1)>=4).astype(int).sum(axis=1)
    counts3 = (circrnas.drop(['circ','genename'], axis=1)>=3).astype(int).sum(axis=1)

    circrnas['counts3']=counts3
    circrnas['counts4']=counts4
    circrnas.to_csv(output_file, sep="\t", index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Compute the number of samples with at least 4 and 3 CCR')
    parser.add_argument('-i', '--input_file',
                        required=True, help='count_file')
    parser.add_argument('-o', '--output_file',
                        required=True, help='output file')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main(args.input_file, args.output_file)
