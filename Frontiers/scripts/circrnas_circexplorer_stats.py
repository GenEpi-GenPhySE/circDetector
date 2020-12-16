#!/usr/bin/env python3
#usage: python3 scripts/stats_mapping.py -i samples.tsv -o mapping_stats_samples.tsv

# imports:
import os
import argparse
import pandas as pd
import numpy as np
import re
import sys
from natsort import natsorted

from collections import defaultdict


# Utility functions
def create_folder(path):
    if not os.path.exists(path):
        os.mkdir(path)

def eprint(*args, **kwargs):
    print(*args,  file=sys.stderr, **kwargs)


def get_valid_keys():
    keys = ['nb input reads', 'readlen', 'uniquely mapped', 'uniquely mapped %',
            'multiple loci', 'multiple loci %', 'too many loci',
            'too many loci %', 'chimeric reads']
    return keys


def get_circexplorer2_columns():
    """
    chrom	Chromosome
    start	Start of circular RNA
    end	End of circular RNA
    name	Circular RNA/Junction reads
    score	Flag of fusion junction realignment
    strand	+ or - for strand
    thickStart	No meaning
    thickEnd	No meaning
    itemRgb	0,0,0
    exonCount	Number of exons
    exonSizes	Exon sizes
    exonOffsets	Exon offsets
    readNumber	Number of junction reads
    circType	Type of circular RNA
    geneName	Name of gene
    isoformName	Name of isoform
    index	Index of exon or intron
    flankIntron Left intron/Right intron
    """
    header = ["chrom", "start", "end", "name", "score", "strand",
              "thickStart", "thickEnd", "itemRgb", "exonCount",
              "exonSizes", "exonOffsets", "readNumber", "circType",
              "geneName", "isoformName", "index", "flankIntron"]
    return header



def get_annotations(sample, root_annot_dir):
    """
    Extract circexplorer2 annotation
    """
    annot_dir = os.path.join(root_annot_dir, sample)
    annot_file = os.path.join(annot_dir, "circexplorer2_circ.txt")
    names = get_circexplorer2_columns()
    annot = pd.read_csv(annot_file, sep='\t', names=names)
    return annot



def get_circ_key(row):
    """
    Construct an unique key for each circular RNA
    """
    chrom = row["chrom"]
    start = row["start"]
    end = row["end"]
    strand = row["strand"]
    return ":".join(map(str,[chrom, start, end, strand]))


def sample_circrnas(sample, species, root_annot_dir):
    annot = get_annotations(sample, root_annot_dir)

    circrnas = defaultdict()
    for index, row in annot.iterrows():
        key = get_circ_key(row)
        circrnas[key] = { 'ccr': row['readNumber'],
                          'circType': row['circType'],
                          'isoformName': row['isoformName']}
    return circrnas

def species_circrnas(samples_circs):
    """
    Compute circRNAs detected in each sample
    """
    circ_set = set()
    circ_annot = defaultdict(lambda: defaultdict())
    for sample in samples_circs:
        circ_keys = samples_circs[sample].keys()
        circ_set.update(circ_keys)
        # Get circ rnas gene names
        for circ in circ_keys:
            if circ not in circ_annot:
                circ_annot[circ] = { 'genename': samples_circs[sample][circ]['isoformName'],
                                         'circType': samples_circs[sample][circ]['circType']}

    species_circ = []
    sorted_samples = natsorted(samples_circs.keys())
    for circ in circ_set:
        circ_per_sample = [circ, circ_annot[circ]['genename'], circ_annot[circ]['circType']]
        for sample in sorted_samples:
            circrnas = samples_circs[sample]
            if circ in circrnas:
                circ_per_sample.append(circrnas[circ]["ccr"])
            else:
                circ_per_sample.append(0)
        species_circ.append(circ_per_sample)
    columns = ['circ', 'genename', 'circType']
    columns.extend(sorted_samples)
    df = pd.DataFrame(species_circ, columns=columns)
    return df

def main(input_file, detection_dir, output_dir):

    samples = pd.read_csv(input_file, sep='\t', comment="#",
                          dtype=str).set_index("sample")

    stats = defaultdict(lambda: defaultdict())
    for index, row in samples.iterrows():
        species = row["species"]
        eprint("Reading %s circeplorer output file" % (index))
        stats[species][index] = sample_circrnas(index, row["species"],
                                    detection_dir)

    species_circ = defaultdict()
    for species in stats:
        circrnas = species_circrnas(stats[species])
        output_file = "%s/%s_circexplorer_counts.tsv" % (output_dir, species)
        circrnas.to_csv(output_file, sep="\t", index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Sample file')
    parser.add_argument('-i', '--input_file',
                        required=True, help='Sample file')
    # parser.add_argument('-m', '--mapping_dir',
    #                     required=True, help='root dir of the mapping results')
    parser.add_argument('-d', '--detection_dir',
                        required=True, help='root dir of the detection results')
    parser.add_argument('-o', '--output_dir',
                        required=True, help='output dir')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main(args.input_file, args.detection_dir, args.output_dir)
