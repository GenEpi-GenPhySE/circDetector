
All the scripts are in python 3
```
module load system/Python-3.6.3
```


### Constructing the sample files

Construct locally bos_taurus.tsv, ovis_aries.tsv, sus_scrofa.tsv sampes files
(WARNING existing files will be overwritten)
```
python scripts/constructsamplefiles.py -i datasets_oct2020.csv -d /work2/genphyse/genepi/arobic/circRNA
```
Those files are used for the detection pipeline and are therefore copied to the followgin directories
```
/work2/genphyse/genepi/arobic/circRNA/bos_taurus
/work2/genphyse/genepi/arobic/circRNA/ovis_aries
/work2/genphyse/genepi/arobic/circRNA/sus_scrofa
```

### Mapping stats

All sample files are merged to produce a sinlge sample file : samples.tsv

```
python scripts/stats_mapping.py -i samples.tsv -d /work2/genphyse/genepi/arobic/circRNA/mapping -o out.tsv
```
Plots were made using  mapping_stats.ipynb
   
# Annotation stats

```
python scripts/stats_annotation.py -i samples.tsv -m /work2/genphyse/genepi/arobic/circRNA/mapping -d /work2/genphyse/genepi/arobic/circRNA/circrna/detection/circrnas -o annotation_stats.tsv
```
Plots were made using annotation_stats.ipynb

     

