#!/usr/bin/env python3
#usage: python3 scripts/stats_mapping.py -i samples.tsv -o mapping_stats_samples.tsv

# imports:
import os
import argparse
import pandas as pd
import numpy as np
import re
import sys
import pysam

from collections import defaultdict

tab = str.maketrans("ACTG", "TGAC")

def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]


def species_shortname(species):
    first, second = species.split("_")
    return first[0]+second[:2]


def eprint(*args, **kwargs):
    print(*args,  file=sys.stderr, **kwargs)


def main(input_file, genome_fasta, species, output_file):

    short_name = species_shortname(species)
    circrnas = pd.read_csv(input_file, sep='\t')
    fasta = pysam.FastaFile(genome_fasta)

    frag_length = 100

    with open(output_file, "w") as fout:
        for index, row in circrnas.iterrows():
            chrom, start_str, end_str, strand = row['circ'].split(":")
            start, end = map(int, [ start_str, end_str])
            # circexplorer coordinates are 0-based half-open
            left = fasta.fetch(chrom, start , start + frag_length)
            right = fasta.fetch(chrom, end - frag_length , end)
            if strand == "+":
                front = left
                rear = right
            else:
                front = reverse_complement_table(right)
                rear = reverse_complement_table(left)
            seq = rear + front
            seq_id = "%s:%s" % (short_name, row["circ"])
            desc = row['genename']
            l_width= 60
            lines = [seq[i:i+l_width] for i in range(0, len(seq), l_width)]
            fastaseq = ">%s %s\n" %(seq_id, desc)
            fastaseq += "\n".join(lines) + "\n"
            fout.write(fastaseq)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Compute the number of samples with at least 4 and 3 CCR')
    parser.add_argument('-i', '--input_file',
                        required=True, help='circrna tsv file')
    parser.add_argument('-g', '--genome_file',
                        required=True, help='genome fasta file')
    parser.add_argument('-s', '--species',
                        required=True, help='species')
    parser.add_argument('-o', '--output_file',
                        required=True, help='output file')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main(args.input_file, args.genome_file, args.species, args.output_file)
