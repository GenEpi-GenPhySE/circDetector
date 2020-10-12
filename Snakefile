# snakemake -p cow-liver-B004/R1.Chimeric.out.junction cow-liver-B004/R2.Chimeric.out.junction --until all

# imports:
import pandas as pd
from snakemake.utils import report

configfile: "config.yaml"

sample_units = pd.read_table(config["samples"]).set_index("sample_unit", drop=False)
samples = set(sample_units["sample"])
tols=[0,2,5,10]
min_cr=5

workdir: config["wdir"]

# report: config["report"]

def sample_chimeric_function_files(wildcards):
    mapdirs = sample_units[sample_units["sample"] == wildcards.sample]['mapdir']
    chim_files = []
    for mapdir in mapdirs:
        chim_files.append(os.path.join(mapdir, "se", wildcards.read,
                          "Chimeric.out.junction"))
    return chim_files


def get_species(wildcards):
    species = wildcards.sample.split("-")[0]
    if species not in config['annotations']:
        print("Warning the species %s is absent from aconfig file" % species)
        exit(1)
    return config['annotations'][species]


# Wildcard constraints
wildcard_constraints:
    sample="|".join(samples),
    reads="R1|R2"


rule all:
    input:
        expand("{sample}/circ_rnas.bed", sample=samples),
        expand("{sample}/annotation_circRNAs.out", sample=samples),
        expand("{sample}/stats_annotation.tsv", sample=samples),
        expand("{sample}/true_intronic_circRNAs.tsv", sample=samples),
        expand("{sample}/true_exonic_circRNAs.tsv", sample=samples),
        expand("{sample}/true_exonic_comparison.tsv", sample=samples),
        expand("{sample}/subexonic_pleg_circRNAs.tsv", sample=samples),
        expand("{sample}/subexonic_meg_circRNAs.tsv", sample=samples),
        "mapping_stat.tsv",
        "stats_annotation_all.tsv",
        "exonic_comparison_all.tsv",
        "notebook.done",
        "cirRNAcounts_cow_exonic.tsv",
        "cirRNAcounts_pig_exonic.tsv"
      

rule formatcountmatrixexonic:
    input:
        config["samples"]
    output:
        "cirRNAcounts_cow_exonic.tsv",
        "cirRNAcounts_pig_exonic.tsv"
    log:
        notebook = "logs/notebooks/FormatCountMatrixExonic.ipynb"
    notebook:
        "notebooks/FormatCountMatrixExonic.ipynb"
    
        
rule formatcountmatrix:
    input:
        config["samples"]
    output:
        touch("notebook.done")     
    log:
        # optional path to the processed notebook
        notebook = "logs/notebooks/FormatCountMatrix.ipynb"
    notebook:
        "notebooks/FormatCountMatrix.ipynb"


rule mergeexoniccomparison:
    input:
        expand("{sample}/true_exonic_comparison.tsv", sample=samples)
    output:
        "exonic_comparison_all.tsv"
    shell:
         "cat {{cow,pig}}-*/true_exonic_comparison.tsv | cut -f1,3,4 |tail -n+2|sort|uniq > {output}" 


rule mergestatannotation:
    input:
        expand("{sample}/stats_annotation.tsv", sample=samples)
    output:
        "stats_annotation_all.tsv"
    shell:
        "cat {input} >> {output}"


rule statannotation:
    input:
        "{sample}/annotation_circRNAs.out"
    output:
        stats = "{sample}/stats_annotation.tsv",
        intronic = "{sample}/true_intronic_circRNAs.tsv",
        exonic = "{sample}/true_exonic_circRNAs.tsv",
        comp_exonic = "{sample}/true_exonic_comparison.tsv",
        sub_exonic_pleg = "{sample}/subexonic_pleg_circRNAs.tsv",
        sub_exonic_meg = "{sample}/subexonic_meg_circRNAs.tsv"

    log:
        stdout = "logs/{sample}_stat_annotation.o",
        stderr = "logs/{sample}_stat_annotation.e"
    shell:
        "python3 ../scripts/stats_annotation.py -i {input} -o_stats {output.stats}"
        " -oi {output.intronic} -oe {output.exonic} -oce {output.comp_exonic} -osepleg {output.sub_exonic_pleg}"
        " -osemeg {output.sub_exonic_meg}"
        " 1>{log.stdout} 2>{log.stderr}"


rule annotation:
    input:
        circ_detected = "{sample}/circ_rnas.bed",
        annot_exon = get_species
    output:
        "{sample}/annotation_circRNAs.out"
    log:
        stdout = "logs/{sample}_annotation.o",
        stderr = "logs/{sample}_annotation.e"
    shell:
        "python3 ../scripts/circRNA_annotation.py -circ {input.circ_detected}"
        " -annot {input.annot_exon} -fmt bed -o {output} 1>{log.stdout} 2>{log.stderr}"


rule detection:
    input:
        R1="{sample}/R1.Chimeric.out.junction",
        R2="{sample}/R2.Chimeric.out.junction"
    output:
        "{sample}/circ_rnas.bed"
    log:
        stdout = "logs/{sample}_detection.o",
        stderr = "logs/{sample}_detection.e"
    shell:
        "python3 ../scripts/circRNA_detection.py -r1 {input.R1} -r2 {input.R2}"
        " -min_cr {min_cr} -fmt bed -o {output} 1>{log.stdout} 2>{log.stderr}"


rule mergechimeric:
    input:
        sample_chimeric_function_files
    output:
        "{sample}/{read}.Chimeric.out.junction"
    shell:
        "cat {input} > {output}"


rule mergemappingstat:
    input:
        config["samples"]
    output:
        "mapping_stat.tsv"
    shell:
        "python3 ../scripts/stats_mapping.py -i {input} -o {output}"