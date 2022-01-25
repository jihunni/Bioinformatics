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

# To convert BigWig format (nucleotide base level score) into ensembl ID format (gene level score)
  Ref: https://www.biostars.org/p/131706/  
  GENCODE [link](https://www.gencodegenes.org/human/)  
  
  To convert bigWig format into bed, 
    Download utility programs from UCSC  
    ```
    http://hgdownload.cse.ucsc.edu/admin/exe/
    ```
    
    

    

    Next, to execute the conversion,
    ```
    $ bigWigToWig data.bigWig data.wig
    $ wig2bed < data.wig > data.bed
    ```
  
      
    Trouble shooting:  
    `Error: SortDetails.cpp, 1006: Unable to create FILE* for temp file: No space left on device. Out of memory.`  
    
    ```
    $ wig2bed --max-mem 10G --sort-tmpdir=/home/jihun/data/hg38_phyloP100way/tmp_sort < hg38.phyloP100way.wig > hg38.phyloP100way.bed
    ```
    
  To convert gtf format into bed format (GENECODE)  
    > RefSeq vs Ensembl vs GENCODE
    > https://bioinformatics.stackexchange.com/questions/21/feature-annotation-refseq-vs-ensembl-vs-gencode-whats-the-difference
    
    ```
    $ wget -qO- ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_21/gencode.v21.annotation.gtf.gz \
      | gunzip -c - \
      | gtf2bed - \
      | grep -w gene - \
      > gencode.v21.genes.bed
    ```
  
  To map phyloP into Ensembl ID,  
    ```
    $ bedmap --echo --echo-map-id-uniq data.bed gencode.v21.genes.bed > data_with_overlapping_gencode_gene_names.bed
    ```
 
