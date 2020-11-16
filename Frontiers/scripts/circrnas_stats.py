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


def get_mapping_stats(sample, species, root_mapping_dir):
    mapping_dir = os.path.join(root_mapping_dir, species, "mapping", sample)
    logfile = os.path.join(mapping_dir, "pe", "Log.final.out")
    stats = parse_log_file(logfile)
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


def get_annotations(sample, root_annot_dir):
    annot_dir = os.path.join(root_annot_dir, sample)
    annot_file = os.path.join(annot_dir, "circexplorer_circ.txt")
    names = get_circexplorer2_columns()
    annot = pd.read_csv(annot_file, sep='\t', names=names)
    return annot


def sample_stats(sample, species, root_annot_dir, root_mapping_dir):
    mapping = get_mapping_stats(sample, species, root_mapping_dir)
    annot = get_annotations(sample, root_annot_dir)

    nb_reads = int(mapping['nb input reads'])
    read_len = int(mapping['readlen'])
    nb_circ = len(annot[annot['readNumber']>=4])
    values = [sample, nb_reads, read_len, nb_circ]

    return values

def get_circ_key(row):
    chrom = row["chrom"]
    start = row["start"]
    end = row["end"]
    strand = row["strand"]
    return ":".join(map(str,[chrom, start, end, strand]))


def sample_circrnas(sample, species, root_annot_dir, root_mapping_dir):
    mapping = get_mapping_stats(sample, species, root_mapping_dir)
    annot = get_annotations(sample, root_annot_dir)

    nb_reads = int(mapping['nb input reads'])
    read_len = int(mapping['readlen'])
    nb_circ = len(annot[annot['readNumber']>=4])
    values = [sample, nb_reads, read_len, nb_circ]

    circrnas = defaultdict()
    for index, row in annot.iterrows():
        if row['circType'] == "circRNA": # only exonic circ RNAs
            key = get_circ_key(row)
            circrnas[key] = { 'ccr': row['readNumber'],
                             'isoformName': row['isoformName']}
    return circrnas


def species_circrnas(samples_circs):
    circ_set = set()
    circ_genenames = defaultdict()
    for sample in samples_circs:
        circ_keys = samples_circs[sample].keys()
        circ_set.update(circ_keys)
        # Get circ rnas gene names
        for circ in circ_keys:
            if circ not in circ_genenames:
                circ_genenames[circ] = samples_circs[sample][circ]['isoformName']

    species_circ = []
    sorted_samples = sorted(samples_circs.keys())
    for circ in circ_set:
        circ_per_sample = [circ, circ_genenames[circ]]
        for sample in sorted_samples:
            circrnas = samples_circs[sample]
            if circ in circrnas:
                circ_per_sample.append(circrnas[circ]["ccr"])
            else:
                circ_per_sample.append(0)
        species_circ.append(circ_per_sample)
    columns = ['circ', 'genename']
    columns.extend(sorted_samples)
    df = pd.DataFrame(species_circ, columns=columns)
    return df

def main(input_file, mapping_dir, detection_dir):

    samples = pd.read_csv(input_file, sep='\t', comment="#",
                          dtype=str).set_index("sample")

    stats = defaultdict(lambda: defaultdict())
    for index, row in samples.iterrows():
        species = row["species"]
        stats[species][index] = sample_circrnas(index, row["species"],
                                    detection_dir, mapping_dir)

    species_circ = defaultdict()
    for species in stats:
        circrnas = species_circrnas(stats[species])
        output_file = "%s_circ_counts.tsv" % species
        circrnas.to_csv(output_file, sep="\t", index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Sample file')
    parser.add_argument('-i', '--input_file',
                        required=True, help='Sample file')
    parser.add_argument('-m', '--mapping_dir',
                        required=True, help='root dir of the mapping results')
    parser.add_argument('-d', '--detection_dir',
                        required=True, help='root dir of the detection results')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main(args.input_file, args.mapping_dir, args.detection_dir)
