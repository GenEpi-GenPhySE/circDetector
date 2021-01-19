# CircRNAs:

**Données :**

Les trois fichiers de comptage :

  - **total RNAseq :** *sus_scrofa_linear_counts.tsv*
  - **mRNAseq :** *raffine.gnexpr.readcount.matrix.tsv*
  - **circRNA :** *sus_scrofa_bsj_gene_counts.tsv* => utilisé simplement pour identifier les gènes produisant des circRNAs

**3 échantillons :**

  - *ssc_testis_1* (FR22JFZ201100831) : animal 31
  - *ssc_testis_2* (FR62ND6201203305) : animal 05
  - *ssc_testis_3* (FR79JCW201106154) : animal 54

**Consignes :**

  - extraire les comptages pour les échantillons mentionnés 
  - normaliser les données (paragraphe normalisation : http://genoweb.toulouse.inra.fr/~faraut/circRNAFrontiers)          
  - faire un plot et une régression linéaire pour chaque animal (nuage de points) : comparer l'expression estimée avec du totalRNAseq à celle estimée avec du mRNAseq
  => Réaliser une régression linéaire entre ces estimations (1 par animal)
  - faire apparaître dans le plot les gènes qui produisent des ARNcirc (fichier *sus_scrofa_bsj_gene_counts.tsv*) 
  - montrer que la relation linéaire n'est pas la même pour ces derniers 
  
**Figures :**
Régression linéaire : comparaison de l'expression du totalRNAseq et du mRNAseq : 

![circrna](https://github.com/ccerutti88/circRNA/blob/master/LevelExpression/circRNA_mRNAseq_vs_totalRNAseq.png)
  
- Abscisse : niveau d'expression en mRNAseq
- Ordonnée : niveau d'espression en total RNAseq
- Points gris : gènes non concernés par la production de circRNAs
- Points marrons : gènes concernés par la production de circRNAs. Le dégradé de marron correspond au junction-ratio.
