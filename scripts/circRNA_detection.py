#!/usr/bin/env python3
# usage: python3 circRNA_detection.py [-h] [-r1 R1_INPUT_FILE] [-r2 R2_INPUT_FILE]
                                    # [-min_cr CR_THRESHOLD] [-tol TOLERANCE]
                                    # [-fmt {bed,gtf}] [-o OUTPUT_FILE]

# Imports
import sys, os, re, argparse, natsort, time
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
from collections import defaultdict
import pandas as pd
import numpy as np
from networkx import (draw, DiGraph, Graph)
import networkx as netwx
import statistics
from statistics import mode
from time import sleep
from tqdm import tqdm
import circRNA as circ


# Utility functions
def eprint(*args, **kwargs):
    print(*args,  file=sys.stderr, **kwargs)
    dsp_control = circ.DisplayControl.instance()
    if dsp_control.verbose:
        print(*args,  file=sys.stderr, **kwargs)


def merge_junctions(df_a):
    """ Merge two dataframes
    :df_r1: dataframe of R1 reads
    :df_r2: dataframe of R2 reads
    :return the merged dataframe
    """
    df_merged = pd.concat(df_a,ignore_index=True)
    return df_merged


def read_junction_file(file, read_type=None):
    """ Read the data, adds the header of the table
        and adds a column with the corresponding read type.
    :file: file to read
    :read_type: read type (R1 or R2)
    :return the dataframe of reads
    """
    names = ["chr_donor", "pos_donor","strand_donor","chr_acceptor",
            "pos_acceptor","strand_acceptor","junction_type",
            "repeat_left","repeat_right","read_name","pos_s1",
            "CIGAR_s1","pos_s2","CIGAR_s2"] # Columns names of the input data
    # Read data as table with pandas:
    df = pd.read_table(file, sep = '\t', names=names, dtype=str,
                       converters = {'pos_acceptor':int, 'pos_donor':int})
    df['read_type'] = read_type # Add column "read" to specify the read type
    return df


def get_valid_circjunctions(df):
    """ Selects valid circular chimeric reads (ccr)
        according to the conditions of the "valid_ccr" function.
    :df: datatframe to filter
    return: the filtered dataframe
    """
    # Parsing dataframe to select valid rows
    # Boolean array with true if the read is a valid ccr:
    valid_ccrs = []
    #for index, row in tqdm(df.iterrows(), total=df.shape[0]):
    for index, row in circ.tqdm_wrapper(df):
        valid_ccrs.append(valid_ccr(row))
    df_circ_cr = df[valid_ccrs]
    ccr_a = [] # Array of circular junctions
    eprint(" "+"\n"+"Generating array of selected ccr..."+"\n")
    #for index, row in tqdm(df_circ_cr.iterrows(), total=df_circ_cr.shape[0]):
    for index, row in circ.tqdm_wrapper(df_circ_cr):
        ccr =  circ.CCR(row) # Create CCR object
        ccr_a.append(ccr) # Add chimeric reads into array
    eprint("\n")
    return ccr_a


def valid_ccr(row):
    """ Select ccr according to several conditions.
    :row: row of the dataframe
    :return boolean "False" for not valid ccr
            or boolean "True" for valid ccr
    """
    # Regular expression of the left dangling CIGAR:
    dangling_left = re.compile('^\d+S\d+M$')
    # Regular expression of the rigth dangling CIGAR:
    dangling_right = re.compile('^\d+M\d+S$')
    if (row['chr_donor']!=row['chr_acceptor']
            or row['strand_donor']!=row['strand_acceptor']):
        return False
    elif abs(row['pos_donor'] - row['pos_acceptor']) > 10000000:
        return False
    else:
        if row['strand_donor'] == '-' :
            if (row['pos_donor']>row['pos_acceptor']
                    or row['pos_s1']>row['pos_s2']):
                return False
            else:
                if (not dangling_left.match(row['CIGAR_s1'])
                        or not dangling_right.match(row['CIGAR_s2'])):
                    return False
                else:
                    return True
        elif row['strand_donor'] == '+' :
            if (row['pos_donor']<row['pos_acceptor']
                    or row['pos_s1']<row['pos_s2']):
                return False
            else:
                if (not dangling_right.match(row['CIGAR_s1'])
                        or not dangling_left.match(row['CIGAR_s2'])):
                    return False
                else:
                    return True


def get_exact_circrnas(circ_junctions, cr_threshold):
    """
    Detects circ RNAs as CR sharing exactly the same circ signature
    chrom start end strand
    """
    eprint("Exact circular RNA detection")
    circ_rna_dict = defaultdict(list)
    for circ_junc in circ_junctions:
        circ_rna_dict[circ_junc.key].append(circ_junc)

    circ_array = []
    for key, cr_array in circ_rna_dict.items():
        if len(cr_array) >= cr_threshold:
            circ_array.append(circ.CircRNA(args.output_file_format, ccr_array=cr_array))
    return circ_array


def same_circ_signature(cr_a, cr_b, tol):
    """ Two-by-two comparison of chimeric reads
    :cr_a: chimeric read a
    :cr_b: chimeric read b
    :tol: number of different start-end positions tolerated
    :return boolean "False" for not valid conditions
    and "True" for valid conditions
    """
    if cr_a.chrom != cr_b.chrom:
        return False
    if abs(cr_a.start - cr_b.start) > tol:
        return False
    if abs(cr_a.end - cr_b.end) > tol:
        return False
    if cr_a.read_type == cr_b.read_type:
        if cr_a.strand != cr_b.strand:
            return False
    else:
        if cr_a.strand == cr_b.strand:
            return False
    return True


def merge_circ_rnas(component):
    if len(component) == 1:
        return component[0]
    elements = []
    for circ_rna in component:
        elements.extend(circ_rna.ccr_array)
    return circ.CircRNA(args.output_file_format, ccr_array=elements)


def connected_components(circ_a, neighbours, cr_threshold):
    """ Generate an undirected graph and
        cluters ccr into connected components
    :ccr_a: array of ccr object
    :neighbours: dictionary of neighbours chimeric reads
    :return an array of circular RNAs
    """
    undirected = netwx.Graph()
    undirected.add_nodes_from(range(len(circ_a)))
    undirected.add_edges_from(neighbours)
    circ_rnas = []  # Circular RNAs array
    # Parsing of connected components composed of chimeric reads:
    for comp in netwx.connected_components(undirected):
        circ_rna = merge_circ_rnas([circ_a[i] for i in comp])
        if circ_rna.nb_ccr >= cr_threshold:
            circ_rnas.append(circ_rna)
    return circ_rnas


def cumulative_cr(circ_a):
    return sum([c.nb_ccr for c in circ_a])


def compute_fuzzy_circ_rnas(ccr_a, cr_threshold, tolerance):
    """ Get the circular RNA corresponding to each
        connected component of chimeric reads
    :ccr_a: array of ccr object
    :return array of circular RNAs
    """
    if cumulative_cr(ccr_a) < cr_threshold:
        return []
    n = len(ccr_a)
    # Neighbour cr share sign the same circular RNA
    neighbours = []
    circ_rnas = []
    for i in range(n):
        for j in range(i+1, n):
            if ccr_a[j].start > ccr_a[i].end:
                break
            # Two-by-two comparison of chimeric reads:
            if same_circ_signature(ccr_a[i], ccr_a[j], tolerance):
                neighbours.append((i, j))
    circ_rnas = connected_components(ccr_a, neighbours, cr_threshold)
    return circ_rnas


def get_independent_intervals(circ_rnas):
    # Get independant intervals
    independant_intervals = []
    circ_interval = []
    chrom = None
    for circ_rna in circ_rnas:
        if not chrom:
            chrom = circ_rna.chrom
            max_end = circ_rna.end
        if circ_rna.chrom != chrom:
            independant_intervals.append(circ_interval)
            circ_interval = []
            chrom = circ_rna.chrom
            max_end = circ_rna.start
        if circ_rna.chrom != chrom or circ_rna.start > max_end:
            independant_intervals.append(circ_interval)
            circ_interval = []
            chrom = circ_rna.chrom
            max_end = circ_rna.start
        max_end = max(max_end, circ_rna.end)
        circ_interval.append(circ_rna)
    independant_intervals.append(circ_interval)
    return independant_intervals


def get_fuzzy_circrnas(circ_rnas, cr_threshold, tolerance):
    eprint("Fuzzy circular RNA detection")
    circ_rnas = natsort.natsorted(circ_rnas, key=lambda e: (e.chrom,
                                                            e.start))
    intervals = get_independent_intervals(circ_rnas)
    circ_rnas = []
    for interval in intervals:
        fuzzy_circ_rnas = compute_fuzzy_circ_rnas(interval, cr_threshold,
                                                  tolerance)
        if fuzzy_circ_rnas:
            circ_rnas.extend(fuzzy_circ_rnas)
    return circ_rnas


def circrna_detection(circ_juntcions, cr_threshold, tolerance):
    """ Detecte circular RNAs
    :df_ccr: dataframe of circular chimeric reads
    :return the array of circular RNAs detected
    """
    # First pass, detect exact circular RNAs (no tolerance)
    circ_rnas = get_exact_circrnas(circ_juntcions, cr_threshold)
    eprint("nb of tol 0 circRNAS %d" % len(circ_rnas))#########
    # if tolerance > 0 compute fuzzy circrnas from previous list
    if tolerance > 0:
        circ_rnas = get_fuzzy_circrnas(circ_rnas, cr_threshold, tolerance)
    else:
        circ_rnas = [c for c in circ_rnas if c.nb_ccr >= cr_threshold]
    return circ_rnas


def write_tab_circ_rnas(circ_rnas, outputfile, fmt):
    """ Write the table of circ RNA retained
    :circ_rnas: array of circular RNAs object
    :outputfile: the output file path
    :return the table of circular RNAs retained
    :start Chimeric.out.junction file: last intron base in
    the reference genome order in the circular transcript
    :end Chimeric.out.junction file: first intron base in the
    reference genome order in the circular transcript
    :start GTF: first exon position
    :end GTF: last exon position
    """
    nb_line = 0 # name bed file
    circ_rnas = natsort.natsorted(circ_rnas, key=lambda e: (e.chrom, e.start))
    eprint(" "+"\n"+"Writing the table of circular RNAs..."+"\n")
    with open(outputfile, "w") as fout:
        eprint("\nConvertion the Chimeric.out.junction format to the %s format...\n" % fmt)
        for circ_rna in circ_rnas:
            # Convert Chimeric.out.junction format to out format:
            circ_rna.start += 1
            circ_rna.end -= 1
            # Convert out format to the chosen format:
            circ_rna.write_annot(nb_line, fout, fmt)
            nb_line += 1


def stats(min_cr, tol, nb_ccr_tot, nb_ccr, nb_circ_rna_detected, circ_rnas):
          nb_tot_cr_circ = sum([c.nb_ccr for c in circ_rnas])
          return "%d\t%d\t%d\t%d\t%d\t%d" % (min_cr, tol, nb_ccr_tot, nb_ccr,
                                             nb_circ_rna_detected, nb_tot_cr_circ)


def stats2(df):
    return df.groupby('nb_ccr')['nb_ccr'] \
           .count().reset_index(name='nb_circ_rnas')


def main(r1_input_file, r2_input_file, cr_threshold, tolerance,
         output_file):

    eprint("\nReading junction files...\n")
    df_junctions_a = []
    if r1_input_file:
        df_junctions_a.append(read_junction_file(r1_input_file, read_type="R1"))
    if r2_input_file:
        df_junctions_a.append(read_junction_file(r2_input_file, read_type="R2"))

    eprint("Merging junction files...\n")
    df = merge_junctions(df_junctions_a)
    eprint("%d chimeric reads \n" % len(df.index))

    # Filtering data:
    eprint("Selecting valid ccr...\n")
    circ_junctions = get_valid_circjunctions(df)
    eprint("%d circular chimeric reads" % len(circ_junctions))

    # Circular RNAs detection:
    circ_rnas = circrna_detection(circ_junctions, cr_threshold, tolerance) #### 'circRNA' object has no attribute 'ccr_array'
    eprint("%d circRNAs detected" % len(circ_rnas))

    # Write the table of circular RNAs retained:
    circ.write_annotation(circ_rnas, output_file, fmt=args.output_file_format,
                          attributes=["nb_ccr", "genomic_size", "left",
                                      "right", "complete", "nb_distinct"])


def parse_arguments():
    parser = argparse.ArgumentParser(description='Identification of chimeric reads and circular RNAs')
    parser.add_argument('-r1', '--r1_input_file',
                        required=False, help='R1 input file name')
    parser.add_argument('-r2', '--r2_input_file',
                        required=False, help='R2 input file name')
    parser.add_argument('-min_cr', '--cr_threshold', type=int,
                        required=False, default=5,
                        help='Minimum number of chimeric reads')
    parser.add_argument('-tol', '--tolerance', type=int, required=False,
                        default=3,
                        help='Tolerance of different positions to start and end')
    parser.add_argument('-fmt', '--output_file_format',
                        required=False, default='gtf', choices=["bed","gtf"], # bed ou gtf
                        help='Output file format : 0-based (bed) or 1-based (gtf)')
    parser.add_argument('-o', '--output_file', required=False, default='circ_rnas_detected.out',
                        help='Output file path') # gft par défaut # attribut : nb ccr, min_cr_distinct
    parser.add_argument('--verbose', help='Print more info', action='store_true')

    args = parser.parse_args()

    dsp_control = circ.DisplayControl.instance()
    dsp_control.verbose = args.verbose

    if not (args.r1_input_file or args.r2_input_file):
        parser.error('At least one input file is requested')
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main(args.r1_input_file, args.r2_input_file,
         args.cr_threshold, args.tolerance,
         args.output_file)