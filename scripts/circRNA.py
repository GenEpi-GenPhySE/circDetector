#!/usr/bin/env python3

# Imports
from collections import defaultdict
import operator, pandas, re
import pybedtools as pybed
import natsort
from tqdm import tqdm
import pandas as pd


from Singleton import Singleton

@Singleton
class DisplayControl():
     def __init__(self):
         self.verbose = False


def tqdm_wrapper(df):
    display_control = DisplayControl.instance()
    if display_control.verbose:
        return tqdm(df.iterrows(), total=df.shape[0])
    else:
        return df.iterrows()


class CircJunction():
    def __init__(self, chrom, start, end, transcription_strand):
        """Initialize parameters of the CCR object"""
        self.chrom = chrom
        self.start = start
        self.end = end
        self.transcript_strand = transcription_strand

    @property
    def key(self):
        return "-".join(map(str,[self.chrom, self.start, self.end,
                                 self.transcript_strand]))


class CCR(CircJunction):
    def __init__(self, record):
        """
        Initialize parameters of the CCR object from junction record of the
        data frame
        """
        chrom = record['chr_donor']
        start, end, transcript_strand = self._setjunction(record)
        super(CCR, self).__init__(chrom, start, end, transcript_strand)
        self.record = record

    @property
    def CIGAR_s1(self):
        return self.record['CIGAR_s1']

    @property
    def CIGAR_s2(self):
        return self.record['CIGAR_s2']

    @property
    def read_type(self):
        return self.record['read_type']

    def _setjunction(self, record):
        strand = record['strand_donor']  # == record['strand_acptor']
        if strand == "+":
            start = record['pos_acceptor']
            end = record['pos_donor']
        else:
            start = record['pos_donor']
            end = record['pos_acceptor']
        if record['read_type'] == "R2":
            transcript_strand = strand
        else:
            transcript_strand = complement(strand)
        return start, end, transcript_strand


class Annotation():

    # initialise class variable
    counter = 0

    def __init__(self, chrom, start, end, strand, feature, attributes, source,
                 name=None):
        Annotation.counter += 1
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.feature = feature
        self.source = source
        self.attributes = attributes
        if name:
            self.name = name
        else:
            self.name = "circRNA_%d" % Annotation.counter
        if isinstance(attributes, dict):
            if "name" not in attributes:
                attributes["name"] = name
            self.attributes_d = attributes
            self.attribues = "; ".join(["%s = %s " %(key, value) for key, value in attributes.items()])
        else:
            self.attributes += "; name=%s" % self.name
            self.attributes_d = dict(re.split(" |=", item) for item in self.attributes.split("; "))
            self.clean_quotes()


    def _to_gtf_record(self, attributes=["name"]):
        s = [
            self.chrom, self.source, self.feature,
            self.start, self.end, self._get_score(),
            self.strand, self._get_frame(),
            self._get_attributes(attributes)
        ]
        return "\t".join(map(str, s))


    def _to_bed_record(self, attributes=["name"]):
        start_b, end_b = (self.start - 1, self.end)
        s = [
            self.chrom, start_b, end_b,
            self.name, self._get_score(),
            self.strand, self.source,
            self.feature, self._get_phase(),
            self._get_attributes(attributes)
        ]
        return "\t".join(map(str, s))

    def topybed(self):
        start_b, end_b = (self.start - 1, self.end)
        record = [self.chrom, start_b, end_b, self.name, ".", self.strand]
        return pybed.create_interval_from_list(record)

    def topybed_ife(self):
        record = [self.chrom, self.start, self.end, self.name, ".", self.strand]
        return pybed.create_interval_from_list(record)

    def _get_attributes(self, attributes):
        attribute_str = "; ".join([ "%s=%s" % (key, self.get_att(key)) for key in attributes]) + ";"
        return attribute_str


    @property
    def biotype(self):
        long_biotype = self.get_att("transcript_biotype")
        biotype_keys = ["protein_coding", "pseudogene", "snRNA", "snoRNA", "miRNA", "lncRNA"]
        biotype_values = ["c", "pseudo", "snc", "sno", "mi", "lnc"]
        succint_biotype = dict(zip(biotype_keys, biotype_values))
        if long_biotype in succint_biotype:
            return succint_biotype[long_biotype]
        else:
            # ribozyme, Y_RNA, scaRNA
            return long_biotype # "misc" = miscRNA: Miscellaneous RNA. A non-coding RNA that cannot be classified.

    @property
    def length(self):
        return self.end - self.start + 1

    def _get_score(self):
        return "."

    def _get_frame(self):
        return "."

    def _get_phase(self):
        return "."

    def get_att(self, key):
        if key not in self.attributes_d:
            print("Warning %s is not a valid attribute" % key)
            exit(1)
        return self.attributes_d[key]

    def clean_quotes(self):
        d = self.attributes_d
        for key in d:
            d[key] = d[key].replace('\"', '')

    def get_annotation_tag(self, end_type):
        return "_".join([self.get_att("exon_id"), end_type, self.biotype])

    def annot_str(self, fmt="gtf", attributes=["name"]):
        if fmt=="gtf":
            return self._to_gtf_record(attributes)
        elif fmt=="bed":
            return self._to_bed_record(attributes)
        else:
            return "ERROR: Unknown format..."

    def __str__(self):
        return self.annot_str()


class CircRNA(Annotation):
    def __init__(self, fmt, *args, **kwargs):
        """ Initialize parameters of the circRNA object
        :ccr_array: array of circular chimeric reads
        """
        for key, value in kwargs.items():###
            if "ccr_array" in kwargs:
                self._init_from_ccr(kwargs.get("ccr_array"))
            if len(args) > 0:
                self._init_default(fmt, args[0], args[1], args[2], args[3], args[4], args[5], name=value) ###
            self.start_annotation = []
            self.end_annotation = []
            self.start_transcript_id = []
            self.end_transcript_id = []
            self.start_gene_id = []
            self.end_gene_id = []
            self.intron_annotation = []
            self.infra_exonic_annotation = []


    def _init_from_ccr(self, ccr_array):
        # Number of distinct CIGARs:
        self.ccr_array = ccr_array  # Array of cicular chimeric reads
        rep_ccr = self._consensus_start_end(ccr_array)
        chrom = rep_ccr.chrom
        start = rep_ccr.start + 1 ###
        end = rep_ccr.end - 1 ###
        strand = rep_ccr.transcript_strand
        feature = "circRNA"
        source = "circRNA_detection"
        attributes = "nb_ccr=%d" % self._compute_nb_cr(ccr_array)
        super(CircRNA, self).__init__(chrom, start, end, strand, feature,
                                      attributes, source) # name=None
        self.source = source
        self.feature = feature

        self.CIGAR_s1 = rep_ccr.CIGAR_s1  # CIGAR of the first segment
        self.CIGAR_s2 = rep_ccr.CIGAR_s2  # CIGAR of the second segment
        self.read_type = rep_ccr.read_type  # Read type (R1 or R2)

        self.attributes_d['nb_ccr'] = self._compute_nb_cr(ccr_array)
        self.attributes_d['nb_distinct'] = compute_distinct_cr(ccr_array)
        self.attributes_d['genomic_size'] = self._compute_genomic_size()
        self.attributes_d['left'] = self._compute_left()
        self.attributes_d['right'] = self._compute_right()
        self.attributes_d['complete'] = self._compute_complete()

    @property
    def nb_ccr(self):
        return self.get_att("nb_ccr")

    @property
    def complete(self):
        return self.get_att("complete")

    @property
    def left(self):
        return self.get_att("left")

    @property
    def right(self):
        return self.get_att("right")

    @property
    def nb_distinct(self):
        return self.get_att("nb_distinct")


    def _init_default(self, chrom, start, end, strand, feature, attributes,
                      source, name):
        super(CircRNA, self).__init__(chrom, start, end, strand,
                                      feature, attributes, source, name)

    @property
    def read_names(self):
        """ Return an array of chimeric reads names """
        reads_str = []
        for cr in self.ccr_array:
            reads_str.append(cr.name)
        return reads_str


    def _compute_nb_cr(self, ccr_array):
        return len(self.ccr_array)

    def _compute_genomic_size(self):
        """ Compute the genomic size """
        return self.end - self.start + 1

    def _compute_left(self):
        """ Compute the coordinates of genomic regions to
            left [chr:start strand] format
        """
        return str(self.chrom)+":"+str(self.start)+self.strand

    def _compute_right(self):
        """ Compute the coordinates of genomic regions to
            right [chr:end strand] format
        """
        return str(self.chrom)+":"+str(self.end)+self.strand

    def _compute_complete(self):
        """ Compute the coordinates of genomic regions to
            complete [chr:start-end strand] format
        """
        return str(self.chrom)+":"+str(self.start)+":"+str(self.end)+self.strand

    def _consensus_start_end(self, ccr_array):
        d_start = defaultdict(list)
        d_end = defaultdict(list)
        for ccr in ccr_array:
            d_start[ccr.start].append(ccr)
            d_end[ccr.end].append(ccr)
        dd_start = defaultdict(list)
        dd_end = defaultdict(list)
        for key_start, value_start in d_start.items():
            dd_start[key_start] = compute_distinct_cr(value_start)
        for key_end, value_end in d_end.items():
            dd_end[key_end] = compute_distinct_cr(value_end)
        max_start = max(dd_start.items(), key=operator.itemgetter(1))[0]
        max_end = max(dd_end.items(), key=operator.itemgetter(1))[0]
        for cr in ccr_array:
            cr.start = max_start
            cr.end = max_end
            return cr

    def add_start_annotation(self, annot):
        self.start_annotation.extend(annot)

    def add_end_annotation(self, annot):
        self.end_annotation.extend(annot)

    def add_start_transcript_id(self, annot):
        self.start_transcript_id.extend(annot)

    def add_end_transcript_id(self, annot):
        self.end_transcript_id.extend(annot)

    def add_start_gene_id(self, annot):
        self.start_gene_id.extend(annot)

    def add_end_gene_id(self, annot):
        self.end_gene_id.extend(annot)

    def get_start_annotation_str(self):
        annotations = []
        for exon_tuple in self.start_annotation:
            annotations.append(exon_tuple[0].get_annotation_tag(exon_tuple[1]))
        return ",".join(annotations)

    def get_end_annotation_str(self):
        annotations = []
        for exon_tuple in self.end_annotation:
            annotations.append(exon_tuple[0].get_annotation_tag(exon_tuple[1]))
        return ",".join(annotations)

    def get_start_transcript_id_str(self):
        annotations = []
        for transcript_tuple in self.start_transcript_id:
            annotations.append(transcript_tuple[0].attributes_d["transcript_id"])
        return ",".join(annotations)

    def get_end_transcript_id_str(self):
        annotations = []
        for transcript_tuple in self.end_transcript_id:
            annotations.append(transcript_tuple[0].attributes_d["transcript_id"])
        return ",".join(annotations)

    def get_start_gene_id_str(self):
        annotations = []
        for gene_tuple in self.start_gene_id:
            annotations.append(gene_tuple[0].attributes_d["gene_id"])
        return ",".join(annotations)

    def get_end_gene_id_str(self):
        annotations = []
        for gene_tuple in self.end_gene_id:
            annotations.append(gene_tuple[0].attributes_d["gene_id"])
        return ",".join(annotations)

    def get_intron_annot_str(self, att):
        annotations = []
        for intron, p_overlap in self.intron_annotation:
            upstream = intron.get_att("up_exon").get_att(att)
            downstream = intron.get_att("down_exon").get_att(att)
            annotations.append("%s-%s" % (upstream, downstream))
        return ",".join(annotations)

  
    def get_intron_annot_gene_str(self, att):
        annotations = []
        for intron, p_overlap in self.intron_annotation:
            upstream = intron.get_att("up_exon").get_att(att)
            strand_up = intron.get_att("up_exon").strand
            biotype_up = intron.get_att("up_exon").biotype
            downstream = intron.get_att("down_exon").get_att(att)
            strand_d = intron.get_att("down_exon").strand
            biotype_d = intron.get_att("down_exon").biotype
            annotations.append("%s_%s_%s-%s_%s_%s" % (upstream, strand_up, biotype_up, 
                                                      downstream, strand_d, biotype_d))
        return ",".join(annotations)


    def annotate_intron(self, intron, p_overlap):
        self.intron_annotation.append((intron, p_overlap))


    def get_infra_exonic_annot_str(self, att):
        annotations = []
        for exon in self.infra_exonic_annotation:
            id_ife = exon.get_att(att)
            annotations.append("%s" % (id_ife))
        return ",".join(annotations)

    
    def get_infra_exonic_gene_annot_str(self, att):
        annotations = []
        for exon in self.infra_exonic_annotation:
            id_ife = exon.get_att(att)
            strand = exon.strand
            biotype = exon.biotype
            annotations.append("%s_%s_%s" % (id_ife, strand, biotype))
        return ",".join(annotations)

    
    def get_infra_exonic_exons_start_end(self):
        annotations = []
        for exon in self.infra_exonic_annotation:
            start_exon_ife = exon.start
            end_exon_ife = exon.end
            annotations.append("%s_%s" % (start_exon_ife, end_exon_ife))
        return ",".join(annotations)


    def annotate_ife(self, exon):
        self.infra_exonic_annotation.append(exon)


def complement(strand):
    """ molecular biology complement """
    return "-" if strand == "+" else "+"


def compute_distinct_cr(ccr_array):
    """ From a list of CR compute the number of distincts CR
    Read type with different alignments and different CIGAR
    :ccr_array: array of circular chimeric reads
    """
    cigar_dict = defaultdict(list)
    for cr in ccr_array:
        key = "-".join(map(str, [cr.read_type, cr.CIGAR_s1]))
        cigar_dict[key].append(cr)
    return len(cigar_dict.keys())


def static_vars(**kwargs):
    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func
    return decorate


@static_vars(counter=0, introns_key_dict=defaultdict())
def get_intron_name(chrom, start, end):
    key = "-".join(map(str, [chrom, start, end]))
    if key not in get_intron_name.introns_key_dict:
        get_intron_name.counter += 1
        name = "INTRON_%d" % get_intron_name.counter
        get_intron_name.introns_key_dict[key] = name
    return get_intron_name.introns_key_dict[key]


@static_vars(introns_name_dict=defaultdict())
def new_intron(upstream_exon, donwstreem_exon):
    chrom = upstream_exon.chrom
    start_i = int(upstream_exon.end) + 1
    end_i = int(donwstreem_exon.start) - 1
    strand = upstream_exon.strand
    name = get_intron_name(chrom, start_i, end_i)
    if name in new_intron.introns_name_dict:
        return None
    new_intron.introns_name_dict[name] = 1
    source = "exon_annot"
    feature = "intron"
    attributes = {"up_exon": upstream_exon, "down_exon":  donwstreem_exon}
    intron = Annotation(chrom, start_i, end_i, strand,
                        feature, attributes, source, name=name)
    return intron


def compute_intronic_positions(d_transcripts):
    """
    returns an array of intron Annotation objects
    """
    introns_a = []
    for transcript_exons in d_transcripts.values():
        # sort by transcript
        sorted_exons = natsort.natsorted(transcript_exons,
                                         key=lambda e: (e.chrom, e.start))
        previous = None
        for exon in sorted_exons:
            if previous:
                intron = new_intron(previous, exon)
                if intron is not None:
                    introns_a.append(intron)
            previous = exon
    return introns_a


def write_annotation(annot_a, outputfile, fmt, attributes):
    with open(outputfile, "w") as fout:
        for record in annot_a:
            fout.write(record.annot_str(fmt, attributes)+"\n")


def read_annotation(file, fmt):
    if fmt=="gtf":
        header = ["chrom", "source", "feature", "start", "end",
                  "score", "strand", "phase", "attributes"]
    elif fmt=="bed":
        header = ["chrom", "start", "end", "gene_id", "score",
                  "strand", "source", "feature", "phase", "attributes"]
    else:
        print("Error: Unknown format...")
    return read_annotation_from_file(file, header)


def read_annotation_from_file(file, header):
    annot_df = pd.read_table(file, sep = '\t', names=header,
                             dtype={"chrom":str, "start":str, "end":str})
    annot_a = []
    for index, row in tqdm_wrapper(annot_df):
        chrom = row.chrom
        start = row.start
        end = row.end
        strand = row.strand
        feature = row.feature
        source = row.source if "source" in header else "."
        attributes = row.attributes
        name = row["gene_id"]
        if feature == "circRNA":
            annot = CircRNA(chrom, start, end, strand, feature, attributes,
                            source, name=name)
        else:
            annot = Annotation(chrom, start, end, strand, feature, attributes,
                               source)
        annot_a.append(annot)
    return annot_a
