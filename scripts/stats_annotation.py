#!/usr/bin/env python3
#usage: python3 scripts/stats_annotation.py -i results_pig_testis_31/annotation_circRNAs_f_0_95.out -o stats_annotation.tsv

# Imports:
import circRNA as circ
import os
import argparse
import pandas as pd
import numpy as np
import re
import sys
import csv


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

def intersection(lst1, lst2): 
    return list(set(lst1) & set(lst2)) 

def get_stats(file, df):
    """ Read the circ_rnas_annotation dataframe and return all statistics to tabular 
        format about circRNAs"""
    nb_tot_exonic = len(df)
    sample = get_sample(file)
    header = list(df.columns.values.tolist())
    
    # Arrays:
    exonic_circ_names = []
    infraexonic_tot_names = []
    monoexonic_circ_names = []
    true_exonic = []
    true_intronic = []
    
    # Counters circRNAs type:
    nb_tot_exonic = 0
    nb_monoexonic = 0
    nb_start_end_exonic = 0
    nb_antisens_exonic = 0
    nb_true_intronic = 0
    nb_infraexonic_tot = 0
    nb_infraexonic_sens = 0
    nb_infraexonic_antisens = 0
    nb_single_annotated_junction = 0

    for index, row in df.iterrows():
        if row.nb_ccr >= 5:
            nb_circ_tot = len(df)
            exon_id_start = row.exons_id_start.split(",")
            exon_start = str(exon_id_start[0])
            exon_start_tmp = exon_start.split("_")
            exon_start = exon_start_tmp[0]
            exon_id_end = row.exons_id_end.split(",")
            exon_end = str(exon_id_end[0])
            exon_end_tmp = exon_end.split("_")
            exon_end = exon_end_tmp[0]
            
            # Exonic circRNAs (True exonics circRNAs):
            if ((len(row.exons_id_start) > 0 or len(row.exons_id_end) > 0)
                and len(row.intron_name) == 0):
                nb_tot_exonic += 1
                exonic_circ_names.append(row.circ_rna_name)
                if row.circ_rna_name in exonic_circ_names:
                    if row.strand == "+": 
                        if ("5" and "+") in row.exons_id_start:
                            if len(row.exons_id_end)==0:
                                nb_single_annotated_junction += 1
                            if ("3" and "+") in row.exons_id_end:
                                nb_start_end_exonic += 1
                                true_exonic.append(row)
                                if ((len(exon_id_start)==1 and len(exon_id_end)==1) and (exon_start == exon_end)):
                                    nb_monoexonic += 1
                                    monoexonic_circ_names.append(row.circ_rna_name)
                        if ("3" and "+") in row.exons_id_end:
                            if len(row.exons_id_start)==0:
                                nb_single_annotated_junction += 1
                        if ("3" and "-") in row.exons_id_start:
                            if ("5" and "-") in row.exons_id_end:
                                nb_antisens_exonic += 1
                    elif row.strand == "-":                         
                        if ("5" and "-") in row.exons_id_end:
                            if len(row.exons_id_start)==0:
                                nb_single_annotated_junction += 1
                            if ("3" and "-") in row.exons_id_start:
                                nb_start_end_exonic += 1
                                true_exonic.append(row)
                                if ((len(exon_id_start)==1 and len(exon_id_end)==1) and (exon_start == exon_end)):
                                    nb_monoexonic += 1
                                    monoexonic_circ_names.append(row.circ_rna_name)
                        if ("3" and "-") in row.exons_id_start:
                            if len(row.exons_id_end)==0:
                                nb_single_annotated_junction += 1
                        if ("5" and "+") in row.exons_id_start:
                            if ("3" and "+") in row.exons_id_end:
                                nb_antisens_exonic += 1

            # Infraexonic circRNAs: 
            if ((len(row.gene_id_ife) > 0) and (row.circ_rna_name not in monoexonic_circ_names)):
                nb_infraexonic_tot += 1
                infraexonic_tot_names.append(row.circ_rna_name)
                if row.circ_rna_name in infraexonic_tot_names:
                    if row.strand == "+":
                        if "+" in row.gene_id_ife:
                            nb_infraexonic_sens += 1
                        elif "-" in row.gene_id_ife:
                            nb_infraexonic_antisens += 1
                    elif row.strand == "-":
                        if "-" in row.gene_id_ife:
                            nb_infraexonic_sens += 1
                        elif "+" in row.gene_id_ife:
                            nb_infraexonic_antisens += 1

            # True intronic circRNAs: 
            if len(row.intron_name) > 0:  
                if row.strand == "+":
                    if (row.end_i - row.end) in range(-5,32):
                        if (row.start - row.start_i) in range(-5,5) or (row.start==row.start_i):
                            nb_true_intronic += 1
                            true_intronic.append(row)
                    elif (row.start == row.start_i and (row.end_i - row.end) > 32): 
                        nb_true_intronic += 1
                        true_intronic.append(row)
                elif row.strand == "-":        
                    if (row.start - row.start_i) in range(-5,32):             
                        if ((row.end - row.end_i) in range(-5,5) or (row.end == row.end_i)): 
                            nb_true_intronic += 1
                            true_intronic.append(row)
                    elif (row.end == row.end_i and (row.start - row.start_i) > 32): 
                        nb_true_intronic += 1
                        true_intronic.append(row)

    # Write exonic table for comparaison between tissues/species:
    stats_exonic = write_comparison_table(sample, true_exonic, "exonic", args.output_comp_exonic_file)
    nb_c = stats_exonic[0]
    nb_nc = stats_exonic[1] 

    # Intronic table to get the proportion of each biotype: 
    write_comparison_table(sample, true_intronic, "intronic", args.output_comp_intronic_file)     

    # Write the true exonic table:
    write_circ_table(true_exonic, header, args.output_exonic_file)

    # Write the true intronic table:
    write_circ_table(true_intronic, header, args.output_intronic_file)

    nb_circ_annotated = nb_start_end_exonic + nb_true_intronic
    nb_circ_non_annotated = nb_circ_tot - (nb_circ_annotated + nb_antisens_exonic + nb_infraexonic_antisens)

    return "\t".join(map(str,[sample, nb_circ_tot, nb_tot_exonic, nb_start_end_exonic, nb_single_annotated_junction, nb_antisens_exonic,
                              nb_monoexonic, nb_infraexonic_tot, nb_infraexonic_sens, nb_infraexonic_antisens,
                              nb_true_intronic, nb_circ_annotated, nb_circ_non_annotated, nb_c, nb_nc]))+"\n"


def write_stat_table(stats, output_file):
    with open(output_file, "w") as fout:
        fout.write(stats)


def write_circ_table(self, header, path, index=None, sep="\t", na_rep='', float_format=None,
                     index_label=None, mode='w', encoding=None, date_format=None, decimal='.'):
    """
    Write a circRNAs list to a tabular-separated file (tsv).
    """
    from pandas.core.frame import DataFrame
    df = DataFrame(self)
    # result is only a string if no path provided, otherwise None
    result = df.to_csv(path, index=index, sep=sep, na_rep=na_rep, float_format=float_format, 
                       header=header, index_label=index_label, mode=mode, encoding=encoding, 
                       date_format=date_format, decimal=decimal)
    if path is None:
        return result


def write_comparison_table(sample, circ_rnas, type, output_file_name):
    # Write exonic table for comparaison between tissues/species:
    header_comp = ["chrom:start-end:strand", "nb_ccr", "gene_id", "biotype"]
    df_circ_rnas = pd.DataFrame(circ_rnas, index=None)
    biotypes = []

    with open(output_file_name, 'w') as fout:
        tsv_writer = csv.writer(fout, delimiter='\t')
        tsv_writer.writerow(header_comp)

        for index, row in df_circ_rnas.iterrows():   
            exons_id_start_a = []
            exons_id_end_a = []

            # Select well-annotated true exonic circRNAs bases on the gene_id:
            # Get gene_id:
            genes_id_start = list(set(list(row.gene_id_start.split(","))))
            genes_id_end = list(set(list(row.gene_id_end.split(","))))
            intersect_gene_id = ", ".join(list(set(genes_id_start).intersection(genes_id_end)))

            if len(intersect_gene_id) > 0:         
                # Get key:
                chrom_start = ":".join(map(str,[row.chrom, row.start]))
                chrom_start_end = "-".join(map(str, [chrom_start , row.end]))
                chrom_start_end_strand = ":".join(map(str, [chrom_start_end , row.strand]))
                # Get the ccr number:
                nb_ccr = row.nb_ccr
                # Get biotype:
                exons_id_start = list(set(list(row.exons_id_start.split(","))))
                exons_id_end = list(set(list(row.exons_id_end.split(","))))
                exons_id_start = list(i.split("_") for i in exons_id_start)
                exons_id_end = list(i.split("_") for i in exons_id_end)
                biotypes_start = list(set(list(i[2] for i in exons_id_start)))
                biotypes_end = list(set(list(i[2] for i in exons_id_end)))
                intersect_biotypes = "".join(list(set(biotypes_start).intersection(biotypes_end)))
                # Write line into the output file:
                s = [chrom_start_end_strand, nb_ccr, intersect_gene_id, intersect_biotypes]
                tsv_writer.writerow(s)
                biotypes.append(intersect_biotypes)

    # Dictionary with key = biotype and value = counter:
    d = {}
    d["sample"] = sample
    for biotype in biotypes:
        d[biotype] = d.get(biotype, 0) + 1    
    values = []
    for key, value in d.items():
        if key != 'c' and key != 'sample':
            values.append(value)
    d["nc"] = sum(values)    

    # df_biotypes = pd.DataFrame(d, index=[0], columns=["sample", "c", "nc"])

    if type=="exonic":
        return d["c"], d["nc"]
    #     df_biotypes = df_biotypes.to_csv(args.output_biotype_file, sep = '\t', index=False)
    

def main():

    # Read the circRNAs annotation file:
    df_circ_annot = read_file(args.input_file)

    # Compute statistics:
    stats = get_stats(args.input_file, df_circ_annot)

    # Write the stats table:
    write_stat_table(stats, args.output_stats_file)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Return annotation statistics tables')
    parser.add_argument('-i', '--input_file', required=True, 
                        help='Annotation circRNAs file')
    parser.add_argument('-o_stats', '--output_stats_file', required=False, 
                        default="stats_annotation.tsv", 
                        help='Table containing statistics of annotation')
    parser.add_argument('-oi', '--output_intronic_file', required=False, 
                        default="true_intronic_circRNAs.tsv",
                        help='Table containing true intronic circRNAs')
    parser.add_argument('-oe', '--output_exonic_file', required=False, 
                        default="true_exonic_circRNAs.tsv",
                        help='Table containing true exonic circRNAs')
    parser.add_argument('-oce', '--output_comp_exonic_file', required=False, 
                        default="true_exonic_comparison.tsv",
                        help='Table containing true exonics circRNAs for comparisons')
    parser.add_argument('-oci', '--output_comp_intronic_file', required=False, 
                        default="true_intronic_comparison.tsv",
                        help='Table containing true exonics circRNAs for comparisons')                     
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_arguments()
    main()