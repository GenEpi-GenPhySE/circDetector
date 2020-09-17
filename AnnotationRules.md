
### The annotation classes

#### Exonic circRNA   (Class-E)
   - Both junctions correspond to exonic boundaries from a single gene.  
   The exonic circRNAs must satisfy the three following rules  
      - the 3' junction of a circRNA must precisely correspond to an exon donnor site (3' end of an exon, ie 5' donnor site of the next intron)  **class-1**    
      - the 5' junction must precisely correspond to an upstream exon acceptor site (5' end of an exon, ie 3' acceptor site of the previous intron)  **class-2**   
      - the exon donor and the exon acceptor are associated to a common gene
      
#### Lariat-derived intronic circRNA and intron circle, intronic circRNA for short  (Class-I)
  - both junctions are located within a single intron:
  - the 5' junction must precisely correspond to the 5' intron donnor site
  - the 3' junction must be compatible with a circularization event limited by the branch point (less than ~60 base pair away from the 3' intron acceptor site) 

#### Sub-exonic circRNA  (Class-SE)
   - both junctions are located within a single exon  **class-3**
   - Currently only the ones that are associated to a gene not reported as lnc, coding gene or pseudo-gene

#### Hypothetic-sub-exonic circRNA    (Class-HSE)
   - The other sub-exonic:
      - both junctions are located within an exon and not of the previous type
 
   [class-1]   The class-1 includes circRNAs with Start-circRNA is identified as known exonic boundary
   
   [class-2]   The class-2 includes circRNAs with End-circRNA is identified as known exonic boundary
   
   [class-3]   The class-3 includes circRNAs mapped inside a single exon 
 
   



### The annotation rules
    For all classes the strand of the circRNA is identical to the strand of the gene

    [class-1]
      for a circRNA located on upstream strand, the Start-circRNA corresponds exactly to an exonic 5'splice site 
           OR   
      for a circRNA located on reverse strand, the Start-circRNA corresponds exactly to an exonic 3'splice site


    [class-2]
     for a circRNA located on upstream strand,the End-circRNA corresponds exactly to an exonic 3' splice site
           OR
     for a circRNA located on reverse strand, the End-circRNA corresponds exactly to an exonic 5' splice site

     Class-E``    
     the circRNA is retained in class-1 and in class-2
     the both associated genes are identical


    [Class-I]
     The circRNA is mapped inside an intron (same strand)
         and
     The 5' boundary of the intron must be compatible with genomic coordinates of the circRNA: -5/5 nt
        and
     The 3' boundary of the intron must be compatible with genomic coordinates of the circRNA: -60/5 nt

    [Class-3]
     The circRNA is mapped inside an exon of a described gene (same strand)
     
    [Class-SE]
     The circRNA is retained in Class-4
         and
     The described gene is not described as 'lnc', 'c', or 'pseudo'

    [Class-HSE]
     The circRNA is retained in Class-4
        and
     The described gene is described as 'lnc', 'c', or 'pseudo'








