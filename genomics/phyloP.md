# PhyloP vs Phastcons
  Ref : https://www.biostars.org/p/1239/
  - PhyloP
    - -log(p-value) under a null hypothesis of neutral evolution.
    - positive value : conservated
    - negative value : faster-than expectation evolution
  - Phastcons
    - a probability that each nucleotide belongs to a conserved element
    - non-linear system. ranged from 0 to 1
    - most of part of genome are 0 score or does not exist. Only small parts have quite high score.
  - the highest conserved regions shows by two score systems sholud be highly correlated.
 
# UCSC Genome Database
  UCSC Genome Database [link](https://hgdownload.soe.ucsc.edu/downloads.html#human)  
  [hg19_phyloP100way](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP100way/)  
  [hg38_phastCons100way](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/)  
  [hg38_phyloP100way](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/)  
    
      
  To download the data
  ```
  $ rsync -avz --progress \
  rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/ ./
  ```
