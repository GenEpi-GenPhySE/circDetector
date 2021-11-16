circDetector (CD) : circular RNAs detection and annotation
===================

# Description:
*circDetector* (CD) is a computational tool for detecting and annotation of circular RNAs (circRNAs) from Total RNA-Seq data. This tool is implemented with the Snakemake workflow management system allowing reproducible and scalable data analyses.
CD was developped to identify circRNAs from reads mapped to the reference genome with the STAR tool (Spliced Transcripts Alignment to a Reference) in Single-End (SE). It consists of identifying reads with a circular junction (*chimeric reads*) from the *chimeric.out.junction* files provided by STAR.
CD also provides a file reporting all statistics of STAR-SE mapping. Mapping informations were extracted from the STAR file *Log.final.out*.

**Note:** The documentation of the CD tool is in progress.

# Table of contents:
1. [Installation](#installation)
2. [Snakemake workflow](#snakemake-worflow)
3. [Example Commands](#example-commands)
4. [Contact & Support](#contact)

## Installation:
The source code can be downloaded from [GitHub](https://github.com/GenEpi-GenPhySE/circRNA.git).

### Compilation and configuration:
To install this tool, you need to first fetch the repository [git repository](https://github.com/ccerutti88/circRNA) or download the corresponding compressed files.

```bash
git clone https://github.com/GenEpi-GenPhySE/circRNA.git
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

### Preparing the sample file:
The script *prepare.py* takes as input an tabulated file containing the paths to all the mapping files and generates a simplified tabular file which will be taken at the entrance of the snakemake workflow.

```bash
python scripts/prepare.py -i metadata.tsv -o samples.tsv
```

* **Input:** *metadata.tsv* is a tabulated file containing the following fields:

|Column|Type  |Description                                                |
|:-----:|:----:|:----------------------------------------------------------:|
|**1**     |string|species (bos_taurus, sus_scrofa)                           |
|2     |string|species_short (bos_taurus = cow, sus_scrofa = pig) 		  |
|3     |string|breed (cow\: Angus, Charolais ; pig\: Yana,Pietrain)         |
|**4**     |string|tissue (testis, liver)  		                              |
|5     |string|sex (male, female)				                          |
|6     |string|age (days, month, years)  					              |
|**7**     |string|animal_name (specific name for each individual)            |
|8     |string|sample_deprecated     				                      |
|9     |string|sample_unit (sample uniq name)                             |
|10    |string|fastq (fastq file name)                                    |
|11    |string|SRA (Sequence Read Archive)                                |
|12    |string|platform (Illumina Hiseq 4000)                             |
|13    |string|technology   						                      |
|**14**    |string|mapdir (file path to mapping files)	                      |

**Note:** The table must contain the following informations (column names): species, tissue, animal_name and mapdir (bold in the table above).

* **Output:** *samples.tsv* is a tabulated file containing the following fields:

|Column|Type  |Description                                                |
|:-----:|:----:|:---------------------------------------------:|
|1     |string|sample (sample name)  |
|2     |string|sample_unit (sample uniq name)								   |
|3     |string|mapdir (file path to mapping files)					   |


## Snakemake worflow:

The Snakefile is composed of the main following rules:
- rule detection : circRNAs detection
- rule mappingstat : analyses of statistics of STAR-SE mapping
- rule annotation : circRNAs annotation

![alt text](https://github.com/ccerutti88/circRNA/blob/master/workflow_CD.png?raw=true)

![alt text](https://github.com/ccerutti88/circRNA/blob/master/workflow_CD.png =100x20)

### Rule detection: circRNAs detection:

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

Example of a command executed by Snakemake:
```bash
python3 scripts/circRNA_detection.py -r1 -r2 -o circ_rnas.bed
```

### Rule mappingstat: circRNAs detection:

CD provides a file reporting all statistics of STAR-SE mapping.
Mapping informations were extracted from the STAR file *Log.final.out*.

* **Input:** Example of a *Log.final.out* file:

```bash
Started job on |       Oct 26 05:21:11
                            Started mapping on |       Oct 26 05:29:19
                                   Finished on |       Oct 26 05:52:39
      Mapping speed, Million of reads per hour |       121.89

                         Number of input reads |       47400458
                     Average input read length |       150
                                   UNIQUE READS:
                  Uniquely mapped reads number |       33941471
                       Uniquely mapped reads % |       71.61%
                         Average mapped length |       149.32
                      Number of splices: Total |       8294615
           Number of splices: Annotated (sjdb) |       7882009
                      Number of splices: GT/AG |       8204363
                      Number of splices: GC/AG |       58632
                      Number of splices: AT/AC |       5362
              Number of splices: Non-canonical |       26258
                     Mismatch rate per base, % |       0.57%
                        Deletion rate per base |       0.02%
                       Deletion average length |       2.11
                       Insertion rate per base |       0.02%
                      Insertion average length |       1.74
                            MULTI-MAPPING READS:
       Number of reads mapped to multiple loci |       3381777
            % of reads mapped to multiple loci |       7.13%
       Number of reads mapped to too many loci |       28860
            % of reads mapped to too many loci |       0.06%
                                 UNMAPPED READS:
 Number of reads unmapped: too many mismatches |       0
      % of reads unmapped: too many mismatches |       0.00%
           Number of reads unmapped: too short |       9943085
                % of reads unmapped: too short |       20.98%
               Number of reads unmapped: other |       18180
                    % of reads unmapped: other |       0.04%
                                 CHIMERIC READS:
                      Number of chimeric reads |       129065
                           % of chimeric reads |       0.27%
```

* **Output:** *reports/mapping_stat.tsv* is a tabulated file containing all the informations of the input file (column) for each sample (row).

### Rule annotation: circRNAs annotation:

To annotate circRNAs, we need to clean short circRNAs manually. Post traitement « manuel » des données pour éliminer les circRNAs trop petits.

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

Example of a command executed by Snakemake:
```bash
python3 scripts/circRNA_annotation.py -circ circ_rnas.bed -annot -o annotation_circRNAs.tsv
```

## Commands Options:

### Example Commands:
Testing the detection rule of the snakemake pipeline for a single sample:

```bash
srun snakemake pig-testis-31/circ_rnas.bed -p --cores 1 &> snake.log
```

**Note:** Make sure the circrnaenv environment is in this folder.

```bash
sbatch circ_rnas.sh
```

## Contact & Support
For any bug report, feature request, or questions please fill out an issue through circRNA_project's [issue page](https://github.com/ccerutti88/circRNA/issues).

## Copyright and License
This software is released under GNU General Public License (v3.0).
