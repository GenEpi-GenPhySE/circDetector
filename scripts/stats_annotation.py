#!/usr/bin/env python3
#usage: python3 scripts/stats_annotation.py -i results_pig_testis_31/annotation_circRNAs_f_0_95.out > stats_annotation.tsv

# Imports:
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

def intersection(lst1, lst2): 
    return list(set(lst1) & set(lst2)) 

def get_stats(file, df):
    """ Read the circ_rnas_annotation dataframe and return all statistics to tabular 
        format about circRNAs"""
    nb_tot_exonic = len(df)
    sample = get_sample(file)
    
    monoexonic_circ_names = []
    exonic_circ_names = []
    intronic_circ_names = []
    lariat_derived_circ_names = []
    tot_intronic_names = []
    possible_intron_derived_circ_names = []
    possible_antisens_intron_derived_names = []
    infraexonic_tot_names = []
    infraexonic_sens_names = []
    infraexonic_antisens_names = []
    infraintronic_tot_names = []
    infraintronic_sens_names = []
    infraintronic_antisens_names = []
    biotype_exonic = []
    
    # Counters circRNAs type:
    nb_tot_exonic = 0
    nb_monoexonic = 0
    nb_true_exonic = 0
    nb_start_end_exonic = 0
    nb_probable_exonic = 0
    nb_antisens_exonic = 0
    nb_true_intronic = 0
    nb_tot_intronic = 0
    nb_infraexonic_tot = 0
    nb_infraexonic_sens = 0
    nb_infraexonic_antisens = 0
    nb_lariat_derived = 0
    nb_intron_circles = 0
    nb_possible_intron_derived = 0
    nb_possible_antisens_intron_derived = 0
    nb_infraintronic_tot = 0
    nb_infraintronic_sens = 0
    nb_infraintronic_antisens = 0

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
            
            # eprint(exon_start_tmp[2])

            # Exonic circRNAs (True exonics circRNAs):
            if ((len(row.exons_id_start) > 0 or len(row.exons_id_end) > 0)
                and len(row.intron_name) == 0):
                nb_tot_exonic += 1
                exonic_circ_names.append(row.circ_rna_name)

                if row.circ_rna_name in exonic_circ_names:
                    
                    if row.strand == "+": 
                        if ("5" and "+") in row.exons_id_start:
                            if len(row.exons_id_end)==0:
                                nb_probable_exonic += 1
                            if ("3" and "+") in row.exons_id_end:
                                nb_start_end_exonic += 1
                                if ((len(exon_id_start)==1 and len(exon_id_end)==1) and (exon_start == exon_end)):
                                    nb_monoexonic += 1
                                    monoexonic_circ_names.append(row.circ_rna_name)
                        if "+" in row.gene_id_i:
                            nb_infraintronic_sens += 1
                            infraintronic_sens_names.append(row.circ_rna_name)
                        if ("3" and "+") in row.exons_id_end:
                            if len(row.exons_id_start)==0:
                                nb_probable_exonic += 1
                        if ("3" and "-") in row.exons_id_start:
                            if ("5" and "-") in row.exons_id_end:
                                nb_antisens_exonic += 1
                        if "-" in row.gene_id_i:
                            nb_infraintronic_sens += 1
                            infraintronic_sens_names.append(row.circ_rna_name)

                    elif row.strand == "-":                         
                        if ("5" and "-") in row.exons_id_end:
                            if len(row.exons_id_start)==0:
                                nb_probable_exonic += 1
                            if ("3" and "-") in row.exons_id_start:
                                nb_start_end_exonic += 1
                                if ((len(exon_id_start)==1 and len(exon_id_end)==1) and (exon_start == exon_end)):
                                    nb_monoexonic += 1
                                    monoexonic_circ_names.append(row.circ_rna_name)
                        if ("3" and "-") in row.exons_id_start:
                            if len(row.exons_id_end)==0:
                                nb_probable_exonic += 1
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
                            infraexonic_sens_names.append(row.circ_rna_name)
                        elif "-" in row.gene_id_ife:
                            nb_infraexonic_antisens += 1
                            infraexonic_antisens_names.append(row.circ_rna_name)
                    elif row.strand == "-":
                        if "-" in row.gene_id_ife:
                            nb_infraexonic_sens += 1
                            infraexonic_sens_names.append(row.circ_rna_name)
                        elif "+" in row.gene_id_ife:
                            nb_infraexonic_antisens += 1
                            infraexonic_antisens_names.append(row.circ_rna_name)


            # Infraintronic circRNAs:
            if (len(row.intron_name)>0 and len(row.exons_id_start)==0 and len(row.exons_id_end)==0 and len(row.gene_id_ife)==0):
                nb_infraintronic_tot += 1
                infraintronic_tot_names.append(row.circ_rna_name)
                if row.circ_rna_name in infraintronic_tot_names:
                    if row.strand == "+":
                        if "+" in row.gene_id_i:
                            nb_infraintronic_sens += 1
                            infraintronic_sens_names.append(row.circ_rna_name)
                        elif "-" in row.gene_id_i:
                            nb_infraintronic_antisens += 1
                            infraintronic_antisens_names.append(row.circ_rna_name)
                    elif row.strand == "-":
                        if "-" in row.gene_id_i:
                            nb_infraintronic_sens += 1
                            infraintronic_sens_names.append(row.circ_rna_name)
                        elif "+" in row.gene_id_i:
                            nb_infraintronic_antisens += 1
                            infraintronic_antisens_names.append(row.circ_rna_name)
        

            # Intronic circRNAs:
            if len(row.intron_name) > 0:  
                # Total intronic circRNAs:
                nb_tot_intronic += 1
                tot_intronic_names.append(row.circ_rna_name)
                # True intronic circRNAs:
                if row.strand == "+":
                    if (row.end_i - row.end) in range(-2,32):
                        if (row.start - row.start_i) in range(-2,2) or (row.start==row.start_i):
                            nb_true_intronic += 1
                            intronic_circ_names.append(row.circ_rna_name)
                    elif (row.start == row.start_i and (row.end_i - row.end) > 32): 
                        nb_true_intronic += 1
                        intronic_circ_names.append(row.circ_rna_name)
                elif row.strand == "-":        
                    if (row.start - row.start_i) in range(-2,32):             
                        if ((row.end - row.end_i) in range(-2,2) or (row.end == row.end_i)):
                            nb_true_intronic += 1
                            intronic_circ_names.append(row.circ_rna_name)
                    elif (row.end == row.end_i and (row.start - row.start_i) > 32): 
                        nb_true_intronic += 1
                        intronic_circ_names.append(row.circ_rna_name)
            
            if row.circ_rna_name in intronic_circ_names:
                if row.strand == "+":
                    if ("3" and "+") in row.exons_id_start:
                        # Possible intron-derived circRNAs:
                        if "+" in row.gene_id_i:
                            if (row.start - row.start_i) in range(-2,2):
                                if len(row.exons_id_end)==0:
                                    nb_possible_intron_derived += 1
                                    possible_intron_derived_circ_names.append(row.circ_rna_name)
                                    # Lariat-derived circRNAs:
                                    if (row.end - row.end_i) in range(5,32):
                                        nb_lariat_derived += 1 
                                        lariat_derived_circ_names.append(row.circ_rna_name)
                                # Intron-circle circRNAs:
                                elif ("5" and "+") in row.exons_id_end:
                                    if (row.end - row.end_i) in range(-2,2):
                                        nb_intron_circles += 1  
                        # Possible antisens intron-derived circRNAs:
                        elif "-" in row.gene_id_i:  
                            if ((("-") in row.exons_id_start) and (len(row.exons_id_end)==0)):
                                if (((row.start - row.start_i) in range(-2,32)) and ((row.end - row.end_i) in range(-2,32))):
                                    nb_possible_antisens_intron_derived += 1
                                    possible_antisens_intron_derived_names.append(row.circ_rna_name)

                elif row.strand == "-":
                    if ("3" and "-") in row.exons_id_end:
                        # Possible intron-derived circRNAs:
                        if "-" in row.gene_id_i:
                            if (row.end - row.end_i) in range(-2,2):
                                if len(row.exons_id_start)==0:
                                    nb_possible_intron_derived += 1
                                    possible_intron_derived_circ_names.append(row.circ_rna_name)   
                                    # Lariat-derived circRNAs: 
                                    if (row.start - row.start_i) in range(5,32):
                                        nb_lariat_derived += 1 
                                        lariat_derived_circ_names.append(row.circ_rna_name)
                                # Intron-circle circRNAs: 
                                elif ("5" and "-") in row.exons_id_start:
                                    if (row.start - row.start_i) in range(-2,2):
                                        nb_intron_circles += 1
                        # Possible antisens intron-derived circRNAs:
                        elif "+" in row.gene_id_i:
                            if ((("+") in row.exons_id_end) and (len(row.exons_id_start)==0)):
                                if (((row.start - row.start_i) in range(-2,32)) and ((row.end - row.end_i) in range(-2,32))):
                                    nb_possible_antisens_intron_derived += 1
                                    possible_antisens_intron_derived_names.append(row.circ_rna_name)
                

    nb_circ_intronic_localization = nb_true_intronic - (nb_possible_intron_derived + nb_intron_circles + nb_possible_antisens_intron_derived)
    nb_true_exonic = nb_tot_exonic - nb_probable_exonic 
    nb_circ_annotated = nb_true_exonic + nb_true_intronic
    nb_circ_classified = nb_probable_exonic + (nb_tot_intronic - nb_true_intronic)
    nb_circ_non_annotated = nb_circ_tot - (nb_circ_annotated + nb_circ_classified)

    eprint("Total circRNAs:", nb_circ_tot)
    eprint("Exonic circRNAs:", nb_tot_exonic, nb_start_end_exonic, nb_probable_exonic, nb_antisens_exonic)
    # eprint("Biotype exonic circRNAs:", len(biotype_exonic))
    eprint("Monoexonic circRNAs:", nb_monoexonic)
    eprint("Infraexonic circRNAs:", nb_infraexonic_tot, nb_infraexonic_sens, nb_infraexonic_antisens)
    eprint("Infraintronic circRNAs:", nb_infraintronic_tot, nb_infraintronic_sens, nb_infraintronic_antisens)
    eprint("Intronic circRNAs:", nb_tot_intronic, nb_true_intronic, nb_possible_intron_derived, nb_lariat_derived, nb_intron_circles, nb_possible_antisens_intron_derived, nb_circ_intronic_localization)
    eprint("Annotated circRNAs:", nb_circ_annotated)
    eprint("Classified circRNAs:", nb_circ_classified)
    eprint("Non_annotated_circRNAs:", nb_circ_non_annotated)
    
    return "\t".join(map(str,[sample, nb_circ_tot, nb_tot_exonic, nb_start_end_exonic, nb_probable_exonic, nb_antisens_exonic, nb_monoexonic, 
                              nb_infraexonic_tot, nb_infraexonic_sens, nb_infraexonic_antisens,
                              nb_infraintronic_tot, nb_infraintronic_sens, nb_infraintronic_antisens,
                              nb_tot_intronic, nb_true_intronic, nb_possible_intron_derived, nb_lariat_derived, nb_intron_circles, nb_possible_antisens_intron_derived, nb_circ_intronic_localization,
                              nb_circ_annotated, nb_circ_classified, nb_circ_non_annotated]))+"\n"

def write_stat_table(stats, output_file):
    with open(output_file, "w") as fout:
        fout.write(stats)

def main():

    # Read the circRNAs annotation file:
    df_circ_annot = read_file(args.input_file)

    # Compute statistics:
    stats = get_stats(args.input_file, df_circ_annot)
    print(stats)

    # Write the stats table:
    write_stat_table(stats, args.output_file)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Sample file')
    parser.add_argument('-i', '--input_file',
                        required=True, help='Sample file')
    parser.add_argument('-o', '--output_file',
                        required=True, help='Sample file')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main()