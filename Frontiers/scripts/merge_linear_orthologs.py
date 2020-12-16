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


def get_gene_counts(count_file, species_id):
    counts = pd.read_csv(count_file, sep='\t')
    counts = counts.rename(columns={'gene_id': species_id})
    return counts.set_index(species_id)

def main(countsA_file, countsB_file, orthologs_file, output):

    orthologs = pd.read_csv(orthologs_file, sep=' ')

    gene_id_A, gene_id_B = orthologs.columns

    counts_A = get_gene_counts(countsA_file, gene_id_A)
    counts_B = get_gene_counts(countsB_file, gene_id_B)

    ortho_A = orthologs.set_index(gene_id_A)
    species_A_ortho = counts_A.join(ortho_A)
    species_A_ortho = species_A_ortho[~species_A_ortho[gene_id_B].isnull()]

    ortho_B =  orthologs.set_index(gene_id_B)
    species_B_ortho = counts_B.join(ortho_B)
    species_B_ortho = species_B_ortho[~species_B_ortho[gene_id_A].isnull()]

    merged = pd.merge(species_A_ortho, species_B_ortho, how='inner', left_index=True, right_on=gene_id_A).drop(columns=[gene_id_A, gene_id_B])
    
    merged.to_csv(output, sep="\t")


def parse_arguments():
    parser = argparse.ArgumentParser(description='Compute the number of samples with at least 4 and 3 CCR')
    parser.add_argument('-a', '--countsA',
                        required=True, help='circrna counts tsv file for genome A')
    parser.add_argument('-b', '--countsB',
                        required=True, help='circrna counts tsv file for genome B')
    parser.add_argument('-c', '--orthologs',
                        required=True, help='orthologs')
    parser.add_argument('-o', '--output_file',
                        required=True, help='output dir')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main(args.countsA, args.countsB, args.orthologs, args.output_file)
