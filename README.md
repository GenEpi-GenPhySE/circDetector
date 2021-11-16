circDetector (CD) : circular RNAs detection and annotation
===================

## Description:
*circDetector* (CD) is a computational tool for detecting and annotation of circular RNAs (circRNAs) from Total RNA-Seq data. This tool is implemented with the Snakemake workflow management system allowing reproducible and scalable data analyses. 
CD was developped to identify circRNAs from reads mapped to the reference genome with the STAR tool (Spliced Transcripts Alignment to a Reference). It consists of identifying reads with a circular junction (*chimeric reads*) from the *chimeric.out.junction* files provided by STAR. 
CD also provides a file reporting all statistics of STAR-SE mapping. Mapping informations were extracted from the file *Log.final.out*.

# Table of contents:
1. [Installation](#installation)
2. [Snakemake workflow](#snakemake-worflow)
3. [Example Commands](#example-commands)
4. [Contact & Support](#contact)

## Installation:
The source code can be downloaded from [GitHub](https://github.com/ccerutti88/circRNA). 

### Compilation and configuration:
To install this tool, you need to first fetch the repository [git repository](https://github.com/ccerutti88/circRNA) or download the corresponding compressed files. 

```bash
git clone --recursive https://github.com/ccerutti88/circRNA.git circRNA
cd circRNA
```

### Setting the environment:
```bash
/tools/python/3.6.3/bin/python3 -m venv circrnaenv
source circrnaenv/bin/activate
pip install Cython
pip install pybedtools
pip install snakemake
pip install pandas
pip install natsort
pip install tqdm
pip install networkx
```

## Snakemake worflow:

The Snakefile is composed of following rules:
- rule mergechimeric : merge *Chimeric.out.junction* files (R1 and R2)
- rule detection : circRNAs detection
- rule mappingstat : analyses of statistics of STAR-SE mapping 
- rule annotation : circRNAs annotation

#### Rule mergechimeric: Preparing the sample file:

The script *prepare.py* takes as input an tabulated file containing the paths to all the mapping files and generates a simplified tabular file which will be taken at the entrance of the snakemake workflow.

```bash
python scripts/prepare.py -i metadata.tsv -o samples.tsv
```

* **Input:** *metadata.py* is a tabulated file containing the following fields:
|Column|Type  |Description                                                |
|-----:|:----:|:----------------------------------------------------------|
|1     |string|Species (bos_taurus, sus_scrofa)                           |
|2     |string|Species_short (bos_taurus = cow, sus_scrofa = pig) 		  |
|3     |string|Breed (cow: Angus, Charolais ; pig: Yana,Pietrain)         |
|4     |string|Tissue (testis, liver)  		                              |
|5     |string|Sex (male, female)				                          |
|6     |string|Age (days, month, years)  					              |
|7     |string|Animal_name (specific name for each individual)            |
|8     |string|Sample_deprecated     				                      |
|9     |string|Sample_unit (sample uniq name)                             |
|10    |string|Fastq (fastq file name)                                    |
|11    |string|SRA (Sequence Read Archive)                                |
|12    |string|Platform (Illumina Hiseq 4000)                             |
|13    |string|Technology   						                      |
|14    |string|Mapdir (file path to mapping files)	                      |

* **Output:** *samples.tsv* is a tabulated file containing the following fields:
|Column|Type  |Description                                                |
|-----:|:----:|:----------------------------------------------------------|
|1     |string|Sample (concatenation of species_short-tissue-animal_name  |
|2     |string|sample_unit 												  |
|3     |string|Mapdir 											          |


#### Rule detection: circRNAs detection:

```bash
  -r1 R1_INPUT_FILE, --r1_input_file R1_INPUT_FILE
                        R1 input file name
  -r2 R2_INPUT_FILE, --r2_input_file R2_INPUT_FILE
                        R2 input file name
  -min_cr CR_THRESHOLD, --cr_threshold CR_THRESHOLD
                        Minimum number of chimeric reads
  -tol TOLERANCE, --tolerance TOLERANCE
                        Tolerance of different positions to start and end
  -fmt {bed,gtf}, --output_file_format {bed,gtf}
                        Output file format : 0-based (bed) or 1-based (gtf)
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output file path
  --verbose             Print more info
```

```bash
python3 scripts/circRNA_detection.py -r1 -r2 -o circ_rnas.bed
```

#### Rule annotation: circRNAs annotation:

```bash
-circ CIRC_RNA_FILE, --circ_rna_file CIRC_RNA_FILE
                        Circular RNA file path
  -fmt {bed,gtf}, --annot_format {bed,gtf}
                        Annotation file format
  -annot ANNOTATION_FILE, --annotation_file ANNOTATION_FILE
                        Annotation file path
  -oe EXONIC_OUTPUT_FILE, --exonic_output_file EXONIC_OUTPUT_FILE
                        Exonic output file path
  -oi INTRONIC_OUTPUT_FILE, --intronic_output_file INTRONIC_OUTPUT_FILE
                        Intronic output file path
  -o OUTPUT, --output OUTPUT
                        Annotation output file path
  --verbose             Print more info
```

```bash
python3 scripts/circRNA_annotation.py -circ circ_rnas.bed -annot -o annotation_circRNAs.tsv
```

## Commands Options:

### Example Commands:
Testing the snakemake pipeline for a single sample:

```bash
srun snakemake pig-testis-31/circ_rnas.bed -p --cores 1 &> snake.log
```

<aside class="notice">
Make sure the circrnaenv environment is in this folder.
</aside>

```bash
sbatch circ_rnas.sh
```

## Contact & Support
For any bug report, feature request, or questions please fill out an issue through circRNA_project's [issue page](https://github.com/ccerutti88/circRNA/issues).

## Copyright and License
This software is released under GNU General Public License (v3.0).
