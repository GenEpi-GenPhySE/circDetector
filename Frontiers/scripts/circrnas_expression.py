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
from gtfparse import read_gtf

from collections import defaultdict

def eprint(*args, **kwargs):
    print(*args,  file=sys.stderr, **kwargs)


def get_count_file(sample):
    #circrnadir = "/work2/genphyse/genepi/arobic/circRNA/circrna/detection/circrnas"
    circrnadir = "/home/tfaraut/Projets/circRNA/FrontiersCircRNAs/analysis/detection/circrnas"
    count_file = "/".join([circrnadir, sample, "ciriquantexplorer/ciricounts.gtf"])
    return count_file


def parse_count_file(sample):
    count_file = get_count_file(sample)

    df = read_gtf(count_file,
                  column_converters={"bsj": float, 'fsj': float, 'junc_ratio': float})
    circrnas = { circ_id: { 'bsj': bsj, 'fsj': fsj, 'gene_id': gene_id, 'junc_ratio': junc_ratio, 'circ_type': circ_type}
                   for (circ_id, gene_id, bsj, fsj, junc_ratio, circ_type)
                   in zip(df['circ_id'], df["gene_id"], df["bsj"], df["fsj"], df["junc_ratio"], df["circ_type"])
                   }
    return circrnas


def species_circrnas(samples_circs, tag='bsj'):
    circ_set = set()
    circ_genes = defaultdict()
    for sample in samples_circs:
        circ_keys = samples_circs[sample].keys()
        circ_set.update(circ_keys)
        # Get circ rnas gene names
        for circ in circ_keys:
            if circ not in circ_genes:
                circ_genes[circ] = { 'gene_id': samples_circs[sample][circ]['gene_id'],
                                     'circ_type': samples_circs[sample][circ]['circ_type']}

    species_circ = []
    sorted_samples = sorted(samples_circs.keys())
    for circ in circ_set:
        circ_per_sample = [circ, circ_genes[circ]['gene_id'], circ_genes[circ]['circ_type']]
        for sample in sorted_samples:
            circrnas = samples_circs[sample]
            if circ in circrnas:
                circ_per_sample.append(circrnas[circ][tag])
            else:
                circ_per_sample.append(0)
        species_circ.append(circ_per_sample)
    columns = [tag, 'gene_id', 'circ_type']
    columns.extend(sorted_samples)
    df = pd.DataFrame(species_circ, columns=columns)
    return df


def get_selected_circ(circrnas, min_ccr=4, min_sample=2):
    circrnas = species_circrnas(circrnas, tag='bsj')
    sample_count = (circrnas.drop(['bsj', 'gene_id', 'circ_type'], axis=1)>=min_ccr).astype(int).sum(axis=1)
    selected = circrnas[sample_count>=min_sample]
    return selected['bsj']


def get_counts(circrnas, tag, min_ccr, min_sample):
    selected_circ = get_selected_circ(circrnas, min_ccr, min_sample)
    counts = species_circrnas(circrnas, tag).set_index(tag)
    selection = counts.loc[selected_circ]
    eprint("%d out of %d" % (len(selection), len(counts)))
    return selection

def main(input_file, output_dir, tag, min_ccr, min_sample):

    if not path.exists(output_dir):
        eprint("Output dir %s is not an existing dir" % output_dir)
        exit(1)

    samples = pd.read_csv(input_file, sep='\t', comment="#").set_index("sample")

    circrnas = defaultdict(lambda: defaultdict())
    for index, row in samples.iterrows():
        species = row["species"]
        eprint("Reading %s count file from %s" % (index, species))
        circrnas[species][index] = parse_count_file(index)

    for species in circrnas:
        eprint("%s" % species)
        counts = get_counts(circrnas[species], tag, min_ccr, min_sample)
        output_file = "%s/%s_%s_counts.tsv" % (output_dir, species, tag)
        counts.to_csv(output_file, sep="\t")


def parse_arguments():
    parser = argparse.ArgumentParser(description='Compute the number of samples with at least 4 and 3 CCR')
    parser.add_argument('-i', '--input_file',
                        required=True, help='circrna tsv file')
    parser.add_argument('-o', '--output_dir',
                        required=True, help='output dir')
    parser.add_argument('-t', '--tag',
                        default='bsj',
                        required=False, help='tag (bsj or fsj, default bsj)')
    parser.add_argument('-c', '--min-ccr',
                        default='4', type=int,
                        required=False, help='tag (bsj or fsj, default bsj)')
    parser.add_argument('-s', '--min-sample',
                        default='2', type=int,
                        required=False, help='tag (bsj or fsj, default bsj)')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main(args.input_file, args.output_dir, args.tag, args.min_ccr, args.min_sample)
