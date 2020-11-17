# Biomart : Get all the orthologous genes between two species
(http://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/)

## Dataset: 
- Ensembl Genes 101
- Pig genes (Sscrofa11.1) / Cow genes (ARS-UCD1.2)

## Filters: 
- **Multi species comparisons:** homologue filters : Orthologous [Sheep/Cow] Genes

## Attributes:
- Homologues
- **Gene:** Gene stable ID
- **Orthologues:**  [Sheep/Cow]  gene stable ID ; [Sheep/Cow] homology type

## Output:
- **head :** gene_ID_[bta/ssc/oar] | gene_ID_[bta/ssc/oar] | homology_type
- *ssc_bta_ortholog_genes.txt*
- *ssc_oar_ortholog_genes.txt*
- *bta_aor_ortholog_genes.txt*

## Select only 1:1 orthologs:

awk -F, '{print $1;$3=="ortholog_one2one"}' ssc_bta_ortholog_genes.txt > ssc_bta_one2one_ortholog_genes.txt
awk '{print $1,$2}' ssc_bta_one2one_ortholog_genes.txt > ssc_bta_one2one_ortholog_genes.tsv
awk -F, '{print $1;$3=="ortholog_one2one"}' ssc_oar_ortholog_genes.txt > ssc_oar_one2one_ortholog_genes.txt
awk '{print $1,$2}' ssc_oar_one2one_ortholog_genes.txt > ssc_oar_one2one_ortholog_genes.tsv
awk -F, '{print $1;$3=="ortholog_one2one"}' bta_aor_ortholog_genes.txt > bta_aor_one2one_ortholog_genes.txt
awk '{print $1,$2}' bta_aor_one2one_ortholog_genes.txt > bta_aor_one2one_ortholog_genes.tsv

## Sort by gene_id:
awk 'NR == 1; NR > 1 {print $0 | "sort -n"}' ssc_bta_one2one_ortholog_genes.tsv > ssc_bta_one2one_ortholog_genes_sort.tsv
awk 'NR == 1; NR > 1 {print $0 | "sort -n"}' ssc_oar_one2one_ortholog_genes.tsv >  ssc_oar_one2one_ortholog_genes_sort.tsv
awk 'NR == 1; NR > 1 {print $0 | "sort -n"}' bta_aor_one2one_ortholog_genes.tsv > bta_aor_one2one_ortholog_genes_sort.tsv
