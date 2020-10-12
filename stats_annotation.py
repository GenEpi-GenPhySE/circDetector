#!/usr/bin/env python3

# Imports:
import os, re, sys, csv, argparse
import circRNA as circ
import pandas as pd
import numpy as np

# Utility functions:
def eprint(*args, **kwargs):
    print(*args,  file=sys.stderr, **kwargs)

def get_sample(file):
    """ Get the sample name from the file path"""
    sample_name = os.path.split(os.path.dirname(file))[1]
    return sample_name

def read_file(file):
    """ Read the circ_rnas_annotation.out file and return a pandas dataframe"""
    df = pd.read_table(file, sep='\t')
    df.replace(np.nan, "", inplace=True)
    return df

def intersection(lst1, lst2): 
    """ Return the intersection between two lists"""
    return list(set(lst1) & set(lst2)) 


def get_exonic_circrnas(df):
    """ Get statistics about exonic circRNAs from the annotation_circRNA.out file"""
    # Arrays, dicts:
    exonic_circ_names, monoexonic_circ_names, true_exonic, single_end_circ_names = [], [], [], []
    ccr_start_end_exonic, ccr_single_annot = dict(), dict()
    # Counter exonic circRNAs:
    nb_monoexonic, nb_start_end_exonic, nb_antisens_exonic, nb_single_annotated_junction = 0, 0, 0, 0

    for index, row in df.iterrows():

        exon_id_start = row.exons_id_start.split(",")
        exon_start = str(exon_id_start[0])
        exon_start_tmp = exon_start.split("_")
        exon_start = exon_start_tmp[0]
        exon_id_end = row.exons_id_end.split(",")
        exon_end = str(exon_id_end[0])
        exon_end_tmp = exon_end.split("_")
        exon_end = exon_end_tmp[0]

        if row.nb_ccr >= 5: 
            if ((len(row.exons_id_start) > 0 or len(row.exons_id_end) > 0)
                and len(row.intron_name) == 0):
                exonic_circ_names.append(row.circ_rna_name)
                if row.circ_rna_name in exonic_circ_names:
                    if row.strand == "+": 
                        if ("5" and "+") in row.exons_id_start:
                            if len(row.exons_id_end)==0:
                                nb_single_annotated_junction += 1
                                ccr_single_annot[row.circ_rna_name] = row.nb_ccr
                                single_end_circ_names.append(row.circ_rna_name)
                            if ("3" and "+") in row.exons_id_end:
                                nb_start_end_exonic += 1
                                true_exonic.append(row)
                                ccr_start_end_exonic[row.circ_rna_name] = row.nb_ccr
                                if ((len(exon_id_start)==1 and len(exon_id_end)==1) and (exon_start == exon_end)):
                                    nb_monoexonic += 1
                                    monoexonic_circ_names.append(row.circ_rna_name)
                        if ("3" and "+") in row.exons_id_end:
                            if len(row.exons_id_start)==0:
                                nb_single_annotated_junction += 1
                                ccr_single_annot[row.circ_rna_name] = row.nb_ccr
                                single_end_circ_names.append(row.circ_rna_name)
                        if ("3" and "-") in row.exons_id_start:
                            if ("5" and "-") in row.exons_id_end:
                                nb_antisens_exonic += 1
                    elif row.strand == "-":                         
                        if ("5" and "-") in row.exons_id_end:
                            if len(row.exons_id_start)==0:
                                nb_single_annotated_junction += 1
                                ccr_single_annot[row.circ_rna_name] = row.nb_ccr
                                single_end_circ_names.append(row.circ_rna_name)
                            if ("3" and "-") in row.exons_id_start:
                                nb_start_end_exonic += 1
                                true_exonic.append(row)
                                ccr_start_end_exonic[row.circ_rna_name] = row.nb_ccr
                                if ((len(exon_id_start)==1 and len(exon_id_end)==1) and (exon_start == exon_end)):
                                    nb_monoexonic += 1
                                    monoexonic_circ_names.append(row.circ_rna_name)
                        if ("3" and "-") in row.exons_id_start:
                            if len(row.exons_id_end)==0:
                                nb_single_annotated_junction += 1
                                ccr_single_annot[row.circ_rna_name] = row.nb_ccr
                                single_end_circ_names.append(row.circ_rna_name)
                        if ("5" and "+") in row.exons_id_start:
                            if ("3" and "+") in row.exons_id_end:
                                nb_antisens_exonic += 1
    return true_exonic, nb_start_end_exonic, nb_single_annotated_junction, single_end_circ_names, ccr_single_annot, monoexonic_circ_names


def get_intronic_circrnas(df):
    # Arrays, dicts:
    intronic_circ_names, true_intronic = [], []
    ccr_true_intronic = dict()
    # Counter intronic circRNAs:
    nb_true_intronic = 0

    for index, row in df.iterrows():
        if row.nb_ccr >= 5: 
            if len(row.intron_name) > 0:  
                if row.strand == "+":
                    if (row.end_i - row.end) in range(-5,60):
                        if (row.start - row.start_i) in range(-5,5) or (row.start==row.start_i):
                            nb_true_intronic += 1
                            true_intronic.append(row)
                            ccr_true_intronic[row.circ_rna_name] = row.nb_ccr
                            intronic_circ_names.append(row.circ_rna_name)
                    elif (row.start == row.start_i and (row.end_i - row.end) > 32): 
                        nb_true_intronic += 1
                        true_intronic.append(row)
                        ccr_true_intronic[row.circ_rna_name] = row.nb_ccr
                        intronic_circ_names.append(row.circ_rna_name)
                elif row.strand == "-":        
                    if (row.start - row.start_i) in range(-5,60):             
                        if ((row.end - row.end_i) in range(-5,5) or (row.end == row.end_i)): 
                            nb_true_intronic += 1
                            true_intronic.append(row)
                            ccr_true_intronic[row.circ_rna_name] = row.nb_ccr
                            intronic_circ_names.append(row.circ_rna_name)
                    elif (row.end == row.end_i and (row.start - row.start_i) > 60): 
                        nb_true_intronic += 1
                        true_intronic.append(row)
                        ccr_true_intronic[row.circ_rna_name] = row.nb_ccr
                        intronic_circ_names.append(row.circ_rna_name)   
    return true_intronic, nb_true_intronic, ccr_true_intronic, intronic_circ_names


def get_subexonic_circrnas(df, monoexonic_circ_names, intronic_circ_names):
    # Subexonic_meg : sens + antisens
    # Subexonic_pleg : sens
    # Arrays, dicts:
    infraexonic_tot_names, infraexonic_circ_names, subexonic, subexonic_antisens_names = [], [], [], []
    ccr_infraexonic = dict()
    # Counter subexonic circRNAs:
    nb_infraexonic, nb_infraexonic_tot, nb_infraexonic_sens, nb_infraexonic_antisens = 0, 0, 0, 0

    for index, row in df.iterrows():
        if row.nb_ccr >= 5: 
            if ((len(row.gene_id_ife) > 0) and (row.circ_rna_name not in monoexonic_circ_names)):
                if (("_5_c_+" not in row.exons_id_start) and ("_5_lnc_+" not in row.exons_id_start) and
                    ("_3_c_-" not in row.exons_id_start) and ("_3_lnc_-" not in row.exons_id_start) and 
                    ("_3_c_+" not in row.exons_id_end) and ("_3_lnc_+" not in row.exons_id_end) and 
                    ("_5_c_-" not in row.exons_id_end) and ("_5_lnc_-" not in row.exons_id_end)): 
                    if len(row.gene_id_start) == 0 and len(row.gene_id_end) == 0:
                        nb_infraexonic_tot += 1
                        infraexonic_tot_names.append(row.circ_rna_name)
                        if row.circ_rna_name not in intronic_circ_names:
                            if row.strand == "+":
                                if "+" in row.gene_id_ife:
                                    nb_infraexonic_sens += 1
                                    subexonic.append(row)
                                    ccr_infraexonic[row.circ_rna_name] = row.nb_ccr
                                    infraexonic_circ_names.append(row.circ_rna_name)
                                elif "-" in row.gene_id_ife:
                                    nb_infraexonic_antisens += 1
                                    subexonic.append(row)
                                    infraexonic_circ_names.append(row.circ_rna_name)
                                    subexonic_antisens_names.append(row.circ_rna_name)
                            elif row.strand == "-":
                                if "-" in row.gene_id_ife:
                                    nb_infraexonic_sens += 1
                                    subexonic.append(row)
                                    infraexonic_circ_names.append(row.circ_rna_name)
                                    ccr_infraexonic[row.circ_rna_name] = row.nb_ccr
                                elif "+" in row.gene_id_ife:
                                    nb_infraexonic_antisens += 1
                                    subexonic.append(row)
                                    infraexonic_circ_names.append(row.circ_rna_name)
                                    subexonic_antisens_names.append(row.circ_rna_name)
                    elif ((len(row.gene_id_start) > 0 and len(row.gene_id_end) == 0) or 
                            (len(row.gene_id_start) == 0  and len(row.gene_id_end) > 0)):
                        nb_infraexonic_tot += 1
                        infraexonic_tot_names.append(row.circ_rna_name)
                        if row.circ_rna_name not in intronic_circ_names:
                            if row.strand == "+":
                                if "+" in row.gene_id_ife:
                                    nb_infraexonic_sens += 1
                                    subexonic.append(row)
                                    infraexonic_circ_names.append(row.circ_rna_name)
                                    ccr_infraexonic[row.circ_rna_name] = row.nb_ccr
                                elif "-" in row.gene_id_ife:
                                    nb_infraexonic_antisens += 1
                                    subexonic.append(row)
                                    infraexonic_circ_names.append(row.circ_rna_name)
                                    subexonic_antisens_names.append(row.circ_rna_name)
                            elif row.strand == "-":
                                if "-" in row.gene_id_ife:
                                    nb_infraexonic_sens += 1
                                    subexonic.append(row)
                                    infraexonic_circ_names.append(row.circ_rna_name)
                                    ccr_infraexonic[row.circ_rna_name] = row.nb_ccr
                                elif "+" in row.gene_id_ife:
                                    nb_infraexonic_antisens += 1
                                    subexonic.append(row)
                                    infraexonic_circ_names.append(row.circ_rna_name)
                                    subexonic_antisens_names.append(row.circ_rna_name)
    return subexonic, infraexonic_circ_names, nb_infraexonic_sens, ccr_infraexonic, subexonic_antisens_names


def get_circrnas(df):
    """ Read the annotation_circRNA.out dataframe and return all statistics to tabular 
        format about circRNAs"""               
    # Get exonic circRNAs:
    stats_exonic_circrnas = get_exonic_circrnas(df)
    exonic = stats_exonic_circrnas[0]
    monoexonic_names = stats_exonic_circrnas[5]
    # Get intronic circRNAs: 
    stats_intronic_circrnas = get_intronic_circrnas(df)
    intronic = stats_intronic_circrnas[0]
    intronic_names = stats_intronic_circrnas[3]
    # Get subexonic circRNAs: 
    stats_subexonic_circrnas = get_subexonic_circrnas(df, monoexonic_names, intronic_names)            
    subexonic = stats_subexonic_circrnas[0]
    subexonic_antisens = stats_subexonic_circrnas[4]
  
    return exonic, subexonic, intronic, subexonic_antisens 


def get_stats_circrnas(sample, df, output_file_name):
    """ Read the annotation_circRNA.out dataframe and return all statistics to tabular 
        format about circRNAs"""
    nb_circ_tot = len(df)  
    # Get circRNAs:
    exonic = get_circrnas(df)[0]    
    subexonic = get_circrnas(df)[1]
    subexonic_antisens = get_circrnas(df)[3]
    intronic = get_circrnas(df)[2]
    
    # Get exonic circRNAs stats:
    stats_exonic = write_comparison_exonic_table(sample, exonic, output_file_name)
    
    # Get subexonic circRNAs stats:
    nb_subexonic = write_subexonic_tables(sample, subexonic, subexonic_antisens, args.output_subexonic_meg_file, 
                                          args.output_subexonic_pleg_file)[0]
    nb_gene_subexonic = write_subexonic_tables(sample, subexonic, subexonic_antisens, args.output_subexonic_meg_file, 
                                               args.output_subexonic_pleg_file)[1]

    if len(stats_exonic) > 5:
        nb_c = stats_exonic[0]
        nb_lnc = stats_exonic[1] 
        nb_autres = stats_exonic[2]
        nb_start_end_exonic_selected = stats_exonic[3]
        nb_ccr_start_end_exonic = stats_exonic[4]
        nb_exonic = write_comparison_exonic_table(sample, exonic, args.output_comp_exonic_file)[5]
    else:
        nb_c = stats_exonic[0] 
        nb_lnc = 0 
        nb_autres = stats_exonic[1]
        nb_start_end_exonic_selected = stats_exonic[2]
        nb_ccr_start_end_exonic = stats_exonic[3]
        nb_exonic = write_comparison_exonic_table(sample, exonic, args.output_comp_exonic_file)[4]

    nb_start_end_exonic = get_exonic_circrnas(df)[1]
    nb_single_annotated_junction = get_exonic_circrnas(df)[2]
    monoexonic_circ_names =  get_exonic_circrnas(df)[5]
    intronic_circ_names = get_intronic_circrnas(df)[3]
    infraexonic_circ_names = get_subexonic_circrnas(df, monoexonic_circ_names, intronic_circ_names)[1]
    single_end_circ_names = get_exonic_circrnas(df)[3]
    nb_infraexonic_sens = get_subexonic_circrnas(df, monoexonic_circ_names, intronic_circ_names)[2]
    nb_true_intronic = get_intronic_circrnas(df)[1]
    ccr_single_annot = get_exonic_circrnas(df)[4]
    ccr_true_intronic = get_intronic_circrnas(df)[2]
    ccr_infraexonic = get_subexonic_circrnas(df, monoexonic_circ_names, intronic_circ_names)[3]

    nb_start_end_false_exonic = nb_start_end_exonic - nb_start_end_exonic_selected
    nb_tot_exonic = nb_start_end_exonic + nb_single_annotated_junction
    nb_common_infra_single = len(intersection(infraexonic_circ_names, single_end_circ_names))
    nb_infraexonic = nb_infraexonic_sens - nb_common_infra_single
    nb_circ_non_annotated = nb_circ_tot - (nb_tot_exonic + nb_infraexonic + nb_true_intronic - nb_start_end_false_exonic)
    nb_ccr_single_annot = sum(ccr_single_annot.values())
    nb_ccr_true_intronic = sum(ccr_true_intronic.values())
    nb_ccr_infraexonic = sum(ccr_infraexonic.values())
    
    return "\t".join(map(str,[sample, nb_circ_tot, nb_exonic, nb_lnc, nb_autres, 
                              nb_subexonic, nb_gene_subexonic, nb_true_intronic]))+"\n"


def write_stats_table(stats, output_file):
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

       
def write_comparison_exonic_table(sample, circ_rnas, output_file_name):
    # Write exonic table for comparaison between tissues/species:
    df_circ_rnas = pd.DataFrame(circ_rnas, index=None)
    header = ["chrom:start-end:strand", "nb_ccr", "gene_id", "biotype", "genomic_size"]
    nb_start_end_annotated, nb_exonic = 0, 0
    nb_ccr_exonics, biotypes = [], []
  
    with open(output_file_name, 'w') as fout:
        tsv_writer = csv.writer(fout, delimiter='\t')
        tsv_writer.writerow(header)
        for index, row in df_circ_rnas.iterrows():  
            exons_id_start_a, exons_id_end_a = [], []
            # Select well-annotated true exonic circRNAs bases on the gene_id:
            # Get key:
            chrom_start = ":".join(map(str,[row.chrom, row.start]))
            chrom_start_end = "-".join(map(str, [chrom_start , row.end]))
            chrom_start_end_strand = ":".join(map(str, [chrom_start_end , row.strand]))
            # Get the ccr number:
            nb_ccr = row.nb_ccr
            nb_ccr_exonics.append(nb_ccr)
            # Get gene_id:
            genes_id_start = list(set(list(row.gene_id_start.split(","))))
            genes_id_end = list(set(list(row.gene_id_end.split(","))))
            # Get biotype:
            exons_id_start = list(set(list(row.exons_id_start.split(","))))
            exons_id_end = list(set(list(row.exons_id_end.split(","))))
            exons_id_start = list(i.split("_") for i in exons_id_start)
            exons_id_end = list(i.split("_") for i in exons_id_end) 
            intersect_gene_id = ", ".join(list(set(genes_id_start).intersection(genes_id_end)))
            if len(intersect_gene_id) > 0: 
                nb_start_end_annotated += 1  
                # Get biotype:    
                biotypes_start = list(set(list(i[2] for i in exons_id_start)))
                biotypes_end = list(set(list(i[2] for i in exons_id_end)))
                intersect_biotypes = "".join(list(set(biotypes_start).intersection(biotypes_end)))
                # Write line into the output file:
                nb_exonic += 1
                genomic_size = row.end - row.start
                s = [chrom_start_end_strand, nb_ccr, intersect_gene_id, intersect_biotypes, genomic_size]
                tsv_writer.writerow(s)
                biotypes.append(intersect_biotypes)  
    # Dictionary with key = biotype and value = counter:
    d = {}
    d["sample"] = sample
    for biotype in biotypes:
        d[biotype] = d.get(biotype, 0) + 1    
    values = []
    valid_keys = ["c", "lnc", "sample"]
    for key, value in d.items():
        if key not in valid_keys:
            values.append(value)
    d["autres"] = sum(values)  
    nb_ccr_exonic = sum(nb_ccr_exonics)
    for key, value in d.items():
        if ("c" in d.keys()) and ("lnc" in d.keys()):
            return d["c"], d["lnc"], d["autres"], nb_start_end_annotated, nb_ccr_exonic, nb_exonic
        elif ("c" in d.keys() and "lnc" not in d.keys()):
            return d["c"], d["autres"], nb_start_end_annotated, nb_ccr_exonic, nb_exonic


def compute_size(exons_start_end):
    pos_start_end = list(exons_start_end.split("_")) 
    pos_start = int(pos_start_end[0])
    pos_end = int(pos_start_end[1])
    size = pos_end - pos_start
    return size 

def get_true_exons_gene_id(exons_start_end, exons_id, gene_id):
    if len(exons_start_end)==1:
        exon_start_end = ",".join(exons_start_end)
        exon_id = ",".join(exons_id)
        gene_id = gene_id
    else:
        annot = dict(zip(exons_start_end, exons_id))
        sizes = []
        for key in annot:
            size = compute_size(key)
            sizes.append(size)
        values = zip(exons_start_end, sizes)
        annot_f = dict(zip(exons_id, values))
        min_size = min([i[1] for i in list(annot_f.values())])
        for key, values in annot_f.items():
            if values[1] == min_size:
                exon_id = key
                exon_start_end = values[0]
    return exon_start_end, exon_id, gene_id


def write_subexonic_tables(sample, circ_rnas, subexonic_antisens, output_file_name_meg, output_file_name_pleg): 
    # Write exonic table for comparaison between tissues/species:
    df_circ_rnas = pd.DataFrame(circ_rnas, index=None)
    nb_sub_exonic = 0
    subexonic_genes = []
    header = list(df_circ_rnas.columns.values.tolist())
  
    with open(output_file_name_meg, 'w') as fout_meg, open(output_file_name_pleg, 'w') as fout_pleg:
        # Subexonic pleg circRNAs:
        tsv_writer_pleg = csv.writer(fout_pleg, delimiter='\t')
        tsv_writer_pleg.writerow(header)
        # Subexonic meg circRNAs:
        tsv_writer_meg = csv.writer(fout_meg, delimiter='\t')
        tsv_writer_meg.writerow(header)

        for index, row in df_circ_rnas.iterrows():  
            # Select well-annotated subexonic circRNAs according to the gene_id:
            # Get gene_id:
            genes_id = list(set(list(row.gene_id_ife.split(","))))
            genes_id = list(i.split("_") for i in genes_id)
            exons_start_end = sorted(set(row.exons_start_end_ife.split(",")), key=row.exons_start_end_ife.split(',').index)
            exon_id = sorted(set(row.exon_id_ife.split(",")), key=row.exon_id_ife.split(',').index)
            gene_id = ",".join(list(set(list(row.gene_id_ife.split(",")))))
            biotypes = []

            exon_start_end = get_true_exons_gene_id(exons_start_end, exon_id, gene_id)[0]
            exon_id = get_true_exons_gene_id(exons_start_end, exon_id, gene_id)[1]
            gene_id = get_true_exons_gene_id(exons_start_end, exon_id, gene_id)[2]

            row.exons_start_end_ife = exon_start_end
            row.exon_id_ife = exon_id
            row.gene_id_ife = gene_id

            if len(genes_id)==1: 
                for i in genes_id:
                    # Subexonic pleg circRNAs:       
                    if ('c' in i or 'lnc' in i or 'pseudo' in i) and (row.circ_rna_name not in subexonic_antisens):
                        nb_sub_exonic += 1
                        s = row
                        tsv_writer_pleg.writerow(s)
                    # Subexonic meg circRNAs:
                    elif 'c' not in i and 'lnc' not in i and 'pseudo' not in i and 'rRNA' not in i:
                        nb_sub_exonic += 1
                        s = row
                        tsv_writer_meg.writerow(s)
                    elif 'rRNA' in i:
                        pass
            
            if len(genes_id) > 1:
                for i in genes_id:
                    if len(i) == 3:
                        biotypes.append(i[2])        
                        biotypes = list(set(biotypes))
                    elif len(i) > 3:
                        p_biotypes = [i[2],i[3]]
                        biotypes = [p for p in p_biotypes]
                if 'rRNA' in biotypes:
                    pass
                # Subexonic pleg circRNAs:
                elif ((biotypes == ['c'] or biotypes == ['lnc'] or biotypes == ['pseudo'] 
                    or ('c' and 'pseudo' in biotypes) or ('c' and 'lnc' in biotypes)
                    or ('lnc' and 'pseudo' in biotypes)) and (row.circ_rna_name not in subexonic_antisens)):
                    nb_sub_exonic += 1
                    s = row
                    tsv_writer_pleg.writerow(s) 
                # Subexonic meg circRNAs:
                else:
                    nb_sub_exonic += 1
                    s = row
                    tsv_writer_meg.writerow(s)

    nb_subexonic_genes = len(list(set([val for sublist in subexonic_genes for val in sublist]))) 
    return nb_sub_exonic, nb_subexonic_genes


def write_circrnas_tables(sample, df, exonic_circrnas, intronic_circrnas, subexonic_circrnas, subexonic_antisens):  
    header = list(df.columns.values.tolist())    
    # Write the exonic circRNAs table:
    write_circ_table(exonic_circrnas, header, args.output_exonic_file)
    # Write the intronic circRNAs table:
    write_circ_table(intronic_circrnas, header, args.output_intronic_file)
    # Write the subexonic circRNAs tables:
    write_subexonic_tables(sample, subexonic_circrnas, subexonic_antisens, args.output_subexonic_meg_file, 
                           args.output_subexonic_pleg_file)


def main():
    # Get sample name:
    sample = get_sample(args.input_file)
    
    # Read the circRNAs annotation file:
    df_circ_annot = read_file(args.input_file)

    # Get the list of exonic, subexonic and intronic circRNAs:
    exonic_circrnas = get_circrnas(df_circ_annot)[0]
    subexonic_circrnas = get_circrnas(df_circ_annot)[1]
    intronic_circrnas = get_circrnas(df_circ_annot)[2]
    subexonic_antisens = get_circrnas(df_circ_annot)[3]

    # Write exonic, intronic and subexonic circRNAs tables:
    circrnas_tables = write_circrnas_tables(sample, df_circ_annot, exonic_circrnas, 
                                            intronic_circrnas, subexonic_circrnas, subexonic_antisens)

    # Compute statistics about exonic, subexonic and intronic circRNAs:
    stats = get_stats_circrnas(sample, df_circ_annot, args.output_comp_exonic_file)

    # Write the circRNAs statistics table:
    write_stats_table(stats, args.output_stats_file)   


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
    parser.add_argument('-osepleg', '--output_subexonic_pleg_file', required=False, 
                        default="subexonic_pleg_circRNAs.tsv",
                        help='Table containing subexonic circRNAs')        
    parser.add_argument('-osemeg', '--output_subexonic_meg_file', required=False, 
                        default="subexonic_meg_circRNAs.tsv",
                        help='Table containing subexonic circRNAs')       
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_arguments()
    main()