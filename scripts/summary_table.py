#!/usr/bin/env python3

# Imports:
import os, re, sys, csv, argparse
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
    """ Read the circRNAs annotation file and return a pandas dataframe"""
    df = pd.read_table(file, sep='\t')
    df.replace(np.nan, "", inplace=True)
    return df


def get_item(item):
    return ''.join(list(set(item)))


def write_summary_meg_table(df, output_file_name, min_size):
    header = ["gene_name", "biotype", "nb_circ_sens", "nb_circ_antisens", "nb_ccr", "exon_name", "start_exon", "end_exon"]
    with open(output_file_name, 'w') as fout:
        tsv_writer = csv.writer(fout, delimiter='\t')
        tsv_writer.writerow(header)
        if len(df) > 0:
            for key, item in df:
                exon_start = get_item(item.exons_start_end_ife).split("_")[0]  
                exon_end = get_item(item.exons_start_end_ife).split("_")[1]
                if int(exon_end) - int(exon_start) > min_size:
                    strand = list(item.strand)
                    gene_id = list(key.split("_"))[0]
                    strand_gene_ife = list(key.split("_"))[1]
                    # Get biotype:
                    if len(list(key.split("_"))) == 3:
                        biotype = list(key.split("_"))[2]
                    else:
                        biotype = list(key.split("_"))[2]+"_"+list(key.split("_"))[3]
                    nb_circ_rna = len(df.get_group(key))
                    nb_circ_sens = strand.count(strand_gene_ife)
                    nb_circ_antisens = nb_circ_rna - nb_circ_sens
                    nb_ccr = sum(item.nb_ccr)
                    exon_name = get_item(item.exon_id_ife)
                    # print(df.get_group(key), "\n\n")
                    row = [gene_id, biotype, nb_circ_sens, nb_circ_antisens, nb_ccr, exon_name, exon_start, exon_end]            
                    tsv_writer.writerow(row)


def write_summary_intronic_table(df, output_file_name, min_size):
    header = ["gene_name", "biotype", "intron_name", "nb_circ_rna", "nb_ccr", "start_intron", "end_intron"]
    with open(output_file_name, 'w') as fout:
        tsv_writer = csv.writer(fout, delimiter='\t')
        tsv_writer.writerow(header)
        if len(df) > 0:
            for key, item in df:
                start_i = ''.join([str(round(x)) for x in list(set(item.start_i))])
                end_i = ''.join([str(round(x)) for x in list(set(item.end_i))])
                if int(end_i) - int(start_i) > min_size:
                    intron_name = key
                    gene_id = get_item(item.gene_id_i).split("_")[0]
                    # Get biotype:
                    biotype = get_item(item.gene_id_i).split(",")
                    if len(biotype)>1:
                        biotype = biotype[0].split("_")[4]
                    else:
                        biotype = ''.join(biotype).split("_")[4]
                    nb_ccr = sum(item.nb_ccr)
                    nb_circ_rna = len(df.get_group(key))
                    row = [gene_id, biotype, intron_name, nb_circ_rna, nb_ccr, start_i, end_i]
                    tsv_writer.writerow(row)


def write_summary_pleg_table(df, output_file_name, min_size):
    header = ["gene_name", "nb_gene", "biotype", "nb_circ_sens", "nb_circ_antisens", "nb_ccr", "exon_name", "start_exon", "end_exon"]
    with open(output_file_name, 'w') as fout:
        tsv_writer = csv.writer(fout, delimiter='\t')
        tsv_writer.writerow(header)
        if len(df) > 0:
            for key, item in df:
                exon_start = get_item(item.exons_start_end_ife).split("_")[0]  
                exon_end = get_item(item.exons_start_end_ife).split("_")[1] 
                if int(exon_end) - int(exon_start) > min_size:
                    # Get gene_id and biotype:
                    gene_id_ife = get_item(item.gene_id_ife).split(",")
                    nb_gene = len(gene_id_ife)
                    if len(gene_id_ife) > 1:
                        gene_id = gene_id_ife[0].split("_")[0]+","+gene_id_ife[1].split("_")[0]
                        strand_gene_ife = get_item(gene_id_ife[0].split("_")[1]+gene_id_ife[1].split("_")[1])
                        biotype = get_item(gene_id_ife[0].split("_")[2]+gene_id_ife[1].split("_")[2])
                    else:
                        gene_id = ''.join(gene_id_ife).split("_")[0]
                        strand_gene_ife = ''.join(gene_id_ife).split("_")[1]
                        biotype = ''.join(gene_id_ife).split("_")[2]
                    strand = list(item.strand)      
                    nb_circ_rna = len(df.get_group(key))
                    nb_circ_sens = strand.count(strand_gene_ife)
                    nb_circ_antisens = nb_circ_rna - nb_circ_sens
                    nb_ccr = sum(item.nb_ccr)
                    exon_name = key
                    row = [gene_id, nb_gene, biotype, nb_circ_sens, nb_circ_antisens, nb_ccr, exon_name, exon_start, exon_end]            
                    tsv_writer.writerow(row)
        

def main():
    
    # Get sample name:
    sample = get_sample(args.input_meg_file)
    min_size = int(args.min_size)

    # Read input annotation circRNAs files:
    df_meg = read_file(args.input_meg_file)
    df_pleg = read_file(args.input_pleg_file)
    df_intronic = read_file(args.input_intronic_file)

    # Meg circRNAs:
    ## Group by gene:
    df_meg = df_meg.groupby('gene_id_ife')
    write_summary_meg_table(df_meg, args.output_meg_file, min_size)

    # Intronic circRNAs:
    ## Group by intron name:
    df_intronic = df_intronic.groupby('intron_name')
    write_summary_intronic_table(df_intronic, args.output_intronic_file, min_size)

    # Pleg circRNAs:
    ## Group by exon:
    df_pleg = df_pleg.groupby('exon_id_ife')
    write_summary_pleg_table(df_pleg, args.output_pleg_file, min_size)

    
def parse_arguments():
    parser = argparse.ArgumentParser(description='Return summary circRNAs tables')
    parser.add_argument('-ip', '--input_pleg_file', required=True, 
                        help='Annotation pleg circRNAs file')
    parser.add_argument('-im', '--input_meg_file', required=True, 
                        help='Annotation meg circRNAs file')
    parser.add_argument('-ii', '--input_intronic_file', required=True, 
                        help='Annotation intronic circRNAs file')
    parser.add_argument('-op', '--output_pleg_file', required=False, 
                        default="pleg_summary.tsv", 
                        help='Table containing summary of pleg circRNAs annotation')     
    parser.add_argument('-om', '--output_meg_file', required=False, 
                        default="meg_summary.tsv", 
                        help='Table containing summary of meg circRNAs annotation')  
    parser.add_argument('-oi', '--output_intronic_file', required=False, 
                        default="intronic_summary.tsv", 
                        help='Table containing summary of intronic circRNAs annotation')
    parser.add_argument('-ms', '--min_size', required=False, default=55, 
                        help='Minimum size of circRNAs')  
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main()