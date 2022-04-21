#!/usr/bin/env python3
#usage: python3 scripts/stats_mapping.py -i samples.tsv -o mapping_stats_samples.tsv

# imports:
import os
import argparse
import pandas as pd
import numpy as np
import re
import sys


# Utility functions
def create_folder(path):
    if not os.path.exists(path):
        os.mkdir(path)

def eprint(*args, **kwargs):
    print(*args,  file=sys.stderr, **kwargs)

def read_file(file):
    """ Read the sample file containing path files and return a pandas dataframe"""
    sample = pd.read_csv(file, sep='\t', dtype=str)
    return sample

def check_mapdirs(mapdirs):
    """ Check if the mapdir exists or not"""
    for mapdir in mapdirs:
        if not os.path.exists(mapdir):
            raise Exception("WARNING following mapdir is missing %s" % mapdir)

def read_log_file(sample):
    """ Read and clean the Log.final.out files containing the summary mapping statistics
        after mapping job is complete and return an array of the pandas dataframe
        containing statistics"""
    stats = []
    d = {}
    reads = ["R1", "R2"]
    paths_out = []
    for index, row in sample.iterrows():
        sample_name = row["sample"]
        sample_unit_name = row["sample_unit"]
        path = row["mapdir"]
        for read in reads:
            path_file = path+"/"+read+"/Log.final.out"
            f = open(path_file, "r")
            for line in f:
                l = line.split("|")
                stats.append(l)
            for stat in stats:
                if len(stat) == 1:
                    stats.remove(stat)
            for stat in stats:
                stat[0] = ' '.join(stat[0].split("\t"))
                stat[0] = ' '.join(stat[0].split())
                stat[1] = ' '.join(stat[1].split())
                stat[1] = stat[1].replace("%", "")
                d['sample'] = sample_name
                d['sample_unit'] = sample_unit_name
                d['read'] = read
                d[stat[0]] = stat[1]
            f.close()
            create_folder("reports/"+sample_unit_name+"/")
            create_folder("reports/"+sample_unit_name+"/"+read+"/")
            path_out = "reports/"+sample_unit_name+"/"+read+"/"+args.output_file
            paths_out.append(path_out)
            write_sample_file(d, path_out)
    return paths_out


def write_sample_file(dict, output_file):
    """Get value from a statistic dictionnary and return a pandas object"""
    df = pd.DataFrame([dict])
    df.to_csv(output_file, sep="\t", index=False)


def write_final_stat_tab(paths_out, output_file):
    df_final = pd.DataFrame()
    for path in paths_out:
        df = pd.read_csv(path, sep="\t")
        df_final = pd.concat([df_final, df], ignore_index=True, sort =False)
    df_final.to_csv(output_file, sep="\t", index=False)


def main():

    # Read the sample file:
    sample = read_file(args.input_file)

    # Check the mapdir of each sample exists or not:
    check_mapdirs(sample["mapdir"])

    create_folder("reports")

    # Read log files and get the path of a table containing statistics:
    paths_out = read_log_file(sample)

    # Write the statistic table with all samples:
    write_final_stat_tab(paths_out, args.output_file)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Sample file')
    parser.add_argument('-i', '--input_file',
                        required=True, help='Sample file')
    parser.add_argument('-o', '--output_file', required=False, 
                        default="mapping_stat.tsv", help='Sample file')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main()
