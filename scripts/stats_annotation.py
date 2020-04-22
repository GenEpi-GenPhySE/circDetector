#!/usr/bin/env python3
#usage: python3 scripts/stats_annotation.py -i results/pig-testis-31/annotation_circRNAs.out 

# imports:
import os
import argparse
import pandas as pd
import numpy as np
import re
import sys


# Utility functions
def eprint(*args, **kwargs):
    print(*args,  file=sys.stderr, **kwargs)

def get_sample(file):
    sample_name = os.path.split(os.path.dirname(file))[1]
    return sample_name

def read_file(file):
    """ Read the circ_rnas_annotation.out file and return a pandas dataframe"""
    df = pd.read_table(file, sep='\t')
    df.replace(np.nan, "", inplace=True)
    return df

def get_stats(file, df):
    """ Read the circ_rnas_annotation dataframe and return all statistics to tabular 
        format about circRNAs"""
    sample = get_sample(file)

    # Counters circRNAs type:
    nb_circ_tot = len(df)
    nb_true_exonic = 0
    nb_probable_exonic = 0
    nb_intronic = 0
    nb_non_annotated = 0
    
    for index, row in df.iterrows():
        if row.nb_ccr >= 5:
            # True exonic circRNAs:
            if (row.strand == "+" and "5" in row.exons_id_start and "3" in row.exons_id_end
                or row.strand == "-" and "3" in row.exons_id_start and "5" in row.exons_id_end):
                nb_true_exonic += 1
            # Probable exonic circRNAs:
            if ((row.strand == "+" and "5" in row.exons_id_start and len(row.exons_id_end)==0) 
                or (row.strand == "+" and "3" in row.exons_id_end and len(row.exons_id_start)==0)
                or (row.strand == "-" and "5" in row.exons_id_end and len(row.exons_id_start)==0)
                or (row.strand == "-" and "3" in row.exons_id_start and len(row.exons_id_end)==0)):
                nb_probable_exonic += 1
            # Intronics circRNAs:
            if len(row.intron_name) != 0:
                nb_intronic += 1
            # Non-annotated:
            if row.exons_id_start == "" and row.exons_id_end == "" and row.intron_name == "":
                nb_non_annotated += 1

    return "\t".join(map(str,[sample, nb_circ_tot, nb_true_exonic, nb_probable_exonic, nb_intronic, nb_non_annotated]))+"\n"


def write_stat_table(stats, output_file):
    with open(output_file, "w") as fout:
        fout.write(stats)


def main():

    # Read the circRNAs annotation file:
    circ_annot = read_file(args.input_file)

    # Compute statistics:
    stats = get_stats(args.input_file, circ_annot)
    
    # Write the stats table:
    write_stat_table(stats, args.output_file)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Sample file')
    parser.add_argument('-i', '--input_file',
                        required=True, help='Sample file')
    parser.add_argument('-o', '--output_file',
                        required=False, help='Sample file')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main()