# flow in genome analysis
![image](https://user-images.githubusercontent.com/48517782/131243597-52de4494-f280-471d-8a0b-ac72d2eed6c9.png)
- data collection : while genome sequencing data, clinical and histological data
- alignment to human genome
  - Genome : hs375d
  - Algorithmn : BWA-Mem
- call somatic variation : subs, indels, SV, CNA
- quality assureance and additional filtering
- downstream scientific analysis

# concepts for genomics
## somatic variants vs germline variants
![image](https://user-images.githubusercontent.com/48517782/131243482-192bc141-fd2e-41c6-b782-265870b9fde0.png)
- germline variants: 
  - Def : variants that are inherited by from the parents via the germ cells, so sperm and oocytes, means the variant has already been present in the genome of at least one of the the parents. 
  - Algorithm: vary against the reference genome (e.g. hg37)
  - Germline variants are either diploid/biallelic, so expected alternative allele frequency is 50% for a heterozygous position.
  ![image](https://user-images.githubusercontent.com/48517782/131243719-9d38bd1a-17ca-4d73-9fdf-8259aa57cf9c.png)  
  Figure. The flowchart of combinations using different sequencers and variant calling pipelines for germline variants.

- somatic variants : 
  - Def : omatic variants arise de novo in the genome of the respective individual. Example: A variant that occurs in a stem cell will be found in all offspring cells that derive from that stem cells, but not in all the other cells of the organism.  
  - Algorithm : (1) different from the control sample (e.g. matched-normal) and (2) different from the reference genome
  - Somatic variants depend on the tumor purity and are not present in all cells tested. As such variant allele frequencies can be much lower.
- In order to distunguish germline from somatic, one sequences the tumor sample and a matched-normal. E.g. in case of lung cancer, one takes the tumor biopsy from the lung, and a matched-normal from the blood. Even though germline variants (risk factor variants) can contribute to pathogenesis, somatic variants are typically more involved a diseases, that is why they are of special interest.
- Without a matched-normal control, one could not distinguish between somatic and germline, because every genome contains tens of thousands of mutations towards the reference genome, so a matched-normal from the same donor is necessary. 
## reference
- https://gatk.broadinstitute.org/hc/en-us/articles/360035890491?id=11127
- https://www.biostars.org/p/305866/

# PCAWG 

![image](https://user-images.githubusercontent.com/48517782/131244392-91fe5eca-0afc-4028-9f16-5bf322c8ef40.png)  
Figure. Scientific output using PCAWG data, in bite-size chunks
