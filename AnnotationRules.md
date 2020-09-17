
### The annotation classes

#### Exonic circRNAs   (Class-E)
   Both junctions correspond to exonic boundaries from a single gene located on the same strand of circRNA
   The exonic circRNAs must satisfy the three following rules  
      - the 3' junction of a circRNA must precisely correspond to an exon donnor site (3' end of an exon, ie 5' donnor site of the next intron)  [class-1]    
      - the 5' junction must precisely correspond to an upstream exon acceptor site (5' end of an exon, ie 3' acceptor site of the previous intron)  [class-2]     
      - the exon donor and the exon acceptor are associated to a common gene
      
#### Single-end Exonic cirRNAs (Class-seE)
   only one junction correspond to an exonic boundary
   
#### Antisens Exonic circRNAs   (Class-asE)
   Both junctions correspond to exonic boundaries from a single gene located on the opposite strand of circRNA  
   The exonic circRNAs must satisfy the three following rules  
      - the 3' junction of a circRNA must precisely correspond to an exon donnor site from a gene located on the opposite strand      
      - the 5' junction must precisely correspond to an upstream exon acceptor site from a gene located on the opposite strand   
      - the exon donor and the exon acceptor are associated to a common gene
      
#### Intronic circRNAs: Lariat-derived intronic circRNA and intron circle  (Class-I)
  - both junctions are located within a single intron:
  - the 5' junction must precisely correspond to the 5' intron donnor site
  - the 3' junction must be compatible with a circularization event limited by the branch point (less than ~60 base pair away from the 3' intron acceptor site) 

#### Sub-exonic circRNAs  (Class-SE)
   - both junctions are located within a single exon  [class-3]
   - Currently only the ones that are associated to a gene not reported as lnc, coding gene or pseudo-gene

#### Hypothetic-sub-exonic circRNAs    (Class-hSE)
   - The other sub-exonic:
      - both junctions are located within an exon and not of the previous type
      
      
      
      
 
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


####     Class-E   
     The circRNA is retained in [class-1] and in [class-2]    
         AND   
     The both associated genes are identical
   
####     Class-seE
     The circRNA is retained in [class-1] OR in [class-2]  BUT NOT in Class-E
     
####     Class-asE   
     For a circRNA located on upstream strand, 
       the Start-circRNA corresponds exactly to an exonic 3'splice site 
         AND
       the End-circRNA corresponds exactly to an exonic 5' splice site
     
     For a circRNA located on reverse strand, 
       the Start-circRNA corresponds exactly to an exonic 5'splice site
         AND
       the End-circRNA corresponds exactly to an exonic 3' splice site

####    Class-I
     The circRNA is mapped inside an intron    
         AND
     the strand of the circRNA is identical to the strand of the gene
         AND
     The 5' boundary of the intron must be compatible with genomic coordinates of the circRNA: -5/5 nt
        AND
     The 3' boundary of the intron must be compatible with genomic coordinates of the circRNA: -60/5 nt

     
 ####   Class-SE
     The circRNA is retained in [class-3]
         AND
     The described gene is not described as 'lnc', 'c', or 'pseudo'


 ####   Class-HSE
     The circRNA is retained in [class-3]
        AND
     The described gene is described as 'lnc', 'c', or 'pseudo'








