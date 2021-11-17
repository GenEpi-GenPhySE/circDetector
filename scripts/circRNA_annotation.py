#!/usr/bin/env python3
# usage: python3 circRNA_annotation.py [-h] [-circ CIRC_RNA_FILE] [-annot ANNOTATION_FILE] 
                                     # [-fmt {bed,gtf}] [-o OUTPUT_FILE]

# Imports
import circRNA as circ
from collections import defaultdict
import sys, os, re, argparse, natsort
import numpy as np
import csv
import pybedtools as pybed

# Utility functions
def eprint(*args, **kwargs):
    print(*args,  file=sys.stderr, **kwargs)
    dsp_control = circ.DisplayControl.instance()
    if dsp_control.verbose:
        print(*args,  file=sys.stderr, **kwargs)


def exonic_annotations(circ_rnas, exons):
    """
    Annotate the circrnas according to the exon annotations
    circ_rnas: an array of cirRNAs object
    exons: an array of circRNAs objects
    returns the circ_rnas array with the circRNAs annotated

    We annotate each circRNA's splice junction with the corresponding
    exons spice junctions using a dict of exons splice junctions
    This means that splice junction correspondances must be exact
    """
    eprint("\nCompute junction dictionary...\n")
    annotation_junctions = compute_junction_dict(exons)
    eprint("\nCircRNA annotation...\n")
    set_annotation(circ_rnas, annotation_junctions)

    return circ_rnas


def junction_key(chrom, pos):
    return "%s:%s" %(chrom, pos)


def compute_junction_dict(annot_array):
    """
    Compute a dictionnary of splic junctions:
    keys: chrom:pos
    values: list of tuples
    tuple: (exon, end_type)
    end_type: 3 (resp 5) for donor (resp acceptor)
    """
    junction_dict = defaultdict(list)
    seen_exons = defaultdict()
    for annot in annot_array:
        if annot.get_att("exon_id") in seen_exons:
            continue
        seen_exons[annot.get_att("exon_id")] = 1
        key_start = junction_key(annot.chrom, annot.start)
        key_end = junction_key(annot.chrom, annot.end)
        if annot.strand == "+":
            start_end_type = "5"
            end_end_type = "3"
        else:
            start_end_type = "3"
            end_end_type = "5"
        # annotating start
        junction_dict[key_start].append((annot, start_end_type))
        # annotating end
        junction_dict[key_end].append((annot, end_end_type))
    return junction_dict


def set_annotation(circ_rnas, annotation_junctions):
    """
    An annotation string would be very much appreciated
    """
    for circ_rna in circ_rnas:
        junction_start = junction_key(circ_rna.chrom, circ_rna.start)
        junction_end = junction_key(circ_rna.chrom, circ_rna.end)
        if junction_start in annotation_junctions:
            circ_rna.add_start_annotation(annotation_junctions[junction_start])
            circ_rna.add_start_transcript_id(annotation_junctions[junction_start])
            circ_rna.add_start_gene_id(annotation_junctions[junction_start])
        if junction_end in annotation_junctions:
            circ_rna.add_end_annotation(annotation_junctions[junction_end])
            circ_rna.add_end_transcript_id(annotation_junctions[junction_end])
            circ_rna.add_end_gene_id(annotation_junctions[junction_end])


def intronic_annotations(circ_rnas, exons):
    """
    Annotate the circrnas according to the introns annotations
    circ_rnas: an array of cirRNAs object
    exons: an array of circRNAs objects
    returns the circ_rnas array with the circRNAs intrno-annotated

    Because at least on circRNA junction is away from existing
    splice junctions, the annotation is solely based on overlap between
    the circRNA and an intron :
      - [start, end] interval must be within the intron
      - this interval must cover at least 90% of the intron
    """

    #eprint("\nCompute introns from exons...\n")
    transcripts = get_exons_per_transcript(exons)
    introns = circ.compute_intronic_positions(transcripts)

    #eprint("\nCompute intersection with circRNAs...\n")
    circ_rnas = annotate_intron_intersection(circ_rnas, introns)

    return circ_rnas


def annotate_intron_intersection(circ_rnas, introns):
    pybed_introns = annot_to_pybed(introns)
    pybed_circrnas = annot_to_pybed(circ_rnas)
    intersections = pybed_circrnas.intersect(pybed_introns, f=0.95, wo=True) 

    circrna_d = {circ.name: circ for circ in circ_rnas}
    introns_d = {intron.name: intron for intron in introns}
    for i in intersections:
        intron = introns_d[i[9]]
        circ = circrna_d[i[3]]
        overlap = i[-1]
        p_overlap = float(overlap)/intron.length # taille circ sur taille intron
        circ.annotate_intron(intron, p_overlap)

    return circ_rnas


def annot_to_pybed(annots):
    pybed_annot = []
    for a in annots:
        pybed_annot.append(a.topybed())
    return pybed.BedTool(pybed_annot).sort()


def infra_exonic_annotations(circ_rnas, exons):
    """
    """
    pybed_exons = annot_to_pybed_ife(exons)
    pybed_circrnas = annot_to_pybed_ife(circ_rnas)  
    intersections = pybed_circrnas.intersect(pybed_exons, f=0.95, wo=True)

    circrna_d = {circ.name: circ for circ in circ_rnas}
    exons_d = {exon.name: exon for exon in exons}
    for i in intersections:
        exon = exons_d[i[9]]
        circ = circrna_d[i[3]]
        circ.annotate_ife(exon)

    return circ_rnas


def annot_to_pybed_ife(annots):
    pybed_annot = []
    for a in annots:
        pybed_annot.append(a.topybed_ife())
    return pybed.BedTool(pybed_annot).sort()


def get_exons_per_transcript(annots):
    """
    returns a dictionnary
    key: ...
    """
    d_transcript = defaultdict(list)
    for exon in annots:
        transcript_id = exon.get_att("transcript_id")
        d_transcript[transcript_id].append(exon)
    return d_transcript


def write_annotated_circrnas(circ_annotated, outputfile):
    circ_rnas = natsort.natsorted(circ_annotated, key=lambda e: (e.chrom, e.start))
    header = "\t".join(["chrom", "start", "end", "strand", "nb_ccr", "circ_rna_name", "exons_id_start", "exons_id_end", 
                        "transcript_id_start", "transcript_id_end", "gene_id_start", "gene_id_end", 
                        "intron_name", "start_i", "end_i", "gene_id_i", "exon_id_i", 
                        'exon_id_ife', "exons_start_end_ife", "gene_id_ife"])
    with open(outputfile, "w") as fout:
        fout.write(header+"\n")
        for circ_rna in circ_rnas:
            intron_annot = circ_rna.intron_annotation 
            if len(intron_annot) == 0:
                start_i = ""
                end_i = ""
                name_i = ""
            else:
                for annot_i in intron_annot:
                    intron = annot_i[0]
                    start_i = intron.start
                    end_i = intron.end
                    name_i = intron.name
            att_d = dict(re.split(" |=", item) for item in circ_rna.attributes.split("; "))
            nb_ccr = att_d["nb_ccr"].replace('\"', '')
            s = [circ_rna.chrom, circ_rna.start, circ_rna.end, 
                circ_rna.strand, nb_ccr, circ_rna.name,
                circ_rna.get_start_annotation_str(), 
                circ_rna.get_end_annotation_str(),
                circ_rna.get_start_transcript_id_str(), 
                circ_rna.get_end_transcript_id_str(), 
                circ_rna.get_start_gene_id_str(), 
                circ_rna.get_end_gene_id_str(),
                name_i, start_i, end_i, 
                circ_rna.get_intron_annot_gene_str('gene_id'),
                circ_rna.get_intron_annot_str('exon_id'), 
                circ_rna.get_infra_exonic_annot_str("exon_id"),
                circ_rna.get_infra_exonic_exons_start_end(),
                circ_rna.get_infra_exonic_gene_annot_str("gene_id")]
            fout.write( "%s\n" % "\t".join(map(str, s)))


def main():

    # Exonic:
    eprint("\nIdentification of EXONIC circRNAs:\n")
    eprint("\nReading circRNA file...\n")
    circ_rnas = circ.read_annotation(args.circ_rna_file, fmt=args.annot_format)
    eprint("\nReading annotation file...\n")
    annots = circ.read_annotation(args.annotation_file, fmt=args.annot_format)

    eprint("\nExonic annotations..\n")
    circ_rnas = exonic_annotations(circ_rnas, annots)

    eprint("\nIntronic annotations..\n")
    circ_rnas = intronic_annotations(circ_rnas, annots)

    #####
    eprint("\nInfra-exonic annotations..\n")
    circ_rnas = infra_exonic_annotations(circ_rnas, annots)
    #####

    eprint("\nWrite annotations..\n")
    write_annotated_circrnas(circ_rnas, args.output)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Annotation')
    parser.add_argument('-circ', '--circ_rna_file',
                        required=True, help='Circular RNA file path')
    parser.add_argument('-fmt', '--annot_format',
                        required=False, default="bed", choices=["bed","gtf"],
                        help='Annotation file format')
    parser.add_argument('-annot', '--annotation_file',
                        required=True, help='Annotataion file path')
    parser.add_argument('-oe', '--exonic_output_file',
                        required=False, default='exonic_circ_rnas_annot.out',
                        help='Exonic output file path')
    parser.add_argument('-oi', '--intronic_output_file',
                        required=False, default='annotation_introns.bed',
                        help='Intronic output file path')
    parser.add_argument('-o', '--output',
                        required=True,
                        help='Annotation output file path')
    parser.add_argument('--verbose', help='Print more info', action='store_true')

    args = parser.parse_args()
    
    dsp_control = circ.DisplayControl.instance()
    dsp_control.verbose = args.verbose
    
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main()
