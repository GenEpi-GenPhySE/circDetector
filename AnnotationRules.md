
### The annotation classes

#### Exonic circRNAs   (Class: Exonic)
   Both junctions correspond to exonic boundaries from a single gene located on the same strand as circRNA
   The circRNAs must satisfy the three following rules  
      - the 3' junction of a circRNA must precisely correspond to an exon donor site (3' end of an exon, ie 5' donor site of the next intron) from a gene located on the same strand as circRNA [class-1]    
      - the 5' junction must precisely correspond to an upstream exon acceptor site (5' end of an exon, ie 3' acceptor site of the previous intron) from a gene located on the same strand as circRNA [class-2]     
      - the exon donor and the exon acceptor are associated to a common gene
      
#### Single-end Exonic cirRNAs (Class: seExo)
   Only one junction of the circRNA corresponds to an exonic boundary   
     or   
  Both junctions correspond to exonic boundaries from two genes
  This class contains no circRNA classified as Exonic
   
#### Antisens Exonic circRNAs   (Class-asExo)
   Both junctions correspond to exonic boundaries from a single gene located on the opposite strand of circRNA  
   The circRNAs must satisfy the three following rules  
      - the 3' junction of a circRNA must precisely correspond to an exon donor site from a gene located on the opposite strand of circRNA     
      - the 5' junction must precisely correspond to an upstream exon acceptor site from a gene located on the opposite strand of circRNA  
      - the exon donor and the exon acceptor are associated to a common gene
    This class contains no circRNA classified as Exonic
      
#### Intronic circRNAs: Lariat-derived intronic circRNA and intron circle  (Class: Intronic)
  - both junctions are located within a single intron
  - the 5' junction must precisely correspond to the 5' intron donor site
  - the 3' junction must be compatible with a circularization event limited by the branch point (less than ~60 base pair away from the 3' intron acceptor site) 

#### Sub-exonic circRNAs from mono-exonic genes  (Class: SubExo-meg)
   - both junctions are located within a single exon  [class-3]
   - only the ones that are associated to a gene not reported as lnc, coding gene or pseudo-gene
      This class contains no circRNA classified as Exonic

#### Sub-exonic circRNAs from pluri-exonic genes   (Class: SubExo-pleg)
   - The other sub-exonic: both junctions are located within an exon 
   - Only the ones that are associated to a gene reported as lnc, coding gene or pseudo-gene
      This class contains no circRNA classified as Exonic
      
      
      
      
 
### The annotation rules
 
 ####    Intermediate classes
    [class-1]  (circRNAs with Start-circRNA is identified as known exonic boundary)
       for a circRNA located on upstream strand, the Start-circRNA corresponds exactly to an exonic 5'splice site 
       for a circRNA located on reverse strand, the Start-circRNA corresponds exactly to an exonic 3'splice site


    [class-2] (circRNAs with End-circRNA is identified as known exonic boundary)
       for a circRNA located on upstream strand,the End-circRNA corresponds exactly to an exonic 3' splice site
       for a circRNA located on reverse strand, the End-circRNA corresponds exactly to an exonic 5' splice site
     
    [Class-3] (circRNAs mapped inside a single exon) 
       The circRNA is mapped inside an exon of a described gene      
         AND  
       The strand of the circRNA is identical to the strand of the gene


####     Class: Exonic   
     The circRNA is retained in [class-1] and in [class-2]    
         AND   
     The both associated genes are identical
   
####     Class: seExo
     The circRNA is retained in [class-1] OR in [class-2]  BUT NOT in Class Exonic
        OR
     The circRNA is retained in [class-1] AND in [class-2]  BUT NOT in Class Exonic
     
####     Class: asExo   
     For a circRNA located on upstream strand, 
       the Start-circRNA corresponds exactly to an exonic 3'splice site 
         AND
       the End-circRNA corresponds exactly to an exonic 5' splice site
         AND
       the both associated genes are identical
     
     For a circRNA located on reverse strand, 
       the Start-circRNA corresponds exactly to an exonic 5'splice site
         AND
       the End-circRNA corresponds exactly to an exonic 3' splice site
         AND
       the both associated genes are identical

####    Class: Intronic
     The circRNA is mapped inside an intron    
         AND
     the strand of the circRNA is identical to the strand of the gene
         AND
     The 5' boundary of the intron must be compatible with genomic coordinates of the circRNA: -5/5 nt
         AND
     The 3' boundary of the intron must be compatible with genomic coordinates of the circRNA: -60/5 nt

     
 ####   Class SubExo-meg
     The circRNA is retained in [class-3]
         AND
     The circRNA is not retained in Class Exonic
         AND
     The described gene is not described as 'lnc', 'c', or 'pseudo'


 ####   Class SubExo-pleg
     The circRNA is retained in [class-3]
         AND
     The circRNA is not retained in Class Exonic
         AND
     The described gene is described as 'lnc', 'c', or 'pseudo'








