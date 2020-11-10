#!/usr/bin/env python3

# imports:
import os
import argparse
import pandas as pd
import numpy as np
import re
import sys
from collections import defaultdict


def read_file(file):
    """ Read the circ_rnas_annotation.out file and return a pandas dataframe"""
    header = ["chrom", "start", "end", "circ_name", "nb_ccr", 
              "strand", "source", "feature", "phase", "attributes"]
    df = pd.read_table(file, sep='\t', names=header)
    df.replace(np.nan, "", inplace=True)
    return df


def get_folder_to_create(df): 
    folders_name = []
    dict_spl_sp = dict(zip(df["sample"], df["species"]))
    split_keys = [x.split("_") for x in dict_spl_sp.keys()]   

    for key in split_keys:
        folders_name.append(key[0]+"_"+key[1])
    folders = list(set(folders_name))
    return folders


def create_folders(folders):
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)


def main():

    # Read samples.tsv file:
    df_samples = pd.read_csv(args.sample_file, sep='\t')

    samples_folders = df_samples["sample"].tolist()
    df_samples = df_samples.loc[df_samples['sample'].isin(samples_folders)] 
    
    # Get the folder list to create:
    folders = get_folder_to_create(df_samples)

    # Create folders:
    create_folders(folders)
   
    # Group folders by species and tissue:
    df_grouped = df_samples.groupby(["species", "tissue"])
    
    for key,item in df_grouped:
        group = df_grouped.get_group(key) 
        # Create col with path of folder:
        group["path_file_bed"] = group["sample"]+"/auzeville.bed"
        
        samples = group["sample"].tolist()
        samples_split = [s.split("_") for s in samples][0]
        folder_name = samples_split[0]+"_"+samples_split[1]
               
        paths = group["path_file_bed"].tolist()

        df_concat = pd.concat(read_file(f) for f in paths)
        df_concat = df_concat.sort_values(by=['chrom', 'start', 'end'])
        df_final = df_concat.groupby(['chrom', 'start', 'end', 'strand'], as_index=False)['nb_ccr'].sum()
        df_final.to_csv(folder_name+"/"+"auzeville.bed", sep="\t", index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Annotation')
    parser.add_argument('-sp', '--sample_file',
                        required=True, help='Sample file path')
    parser.add_argument('-circ', '--circ_rna_file',
                        required=False, help='Circular RNA file path')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_arguments()
    main()