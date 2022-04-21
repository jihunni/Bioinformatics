# reference genome
- Ensembl chromosome other than 1-22, X, and Y | [Biostar](https://www.biostars.org/p/197567/)  
  The KI* and GL* scaffolds are unplaced scaffolds, which often don't have a know chromosome association. In short, these are assembled stretches believed to belong to the human genome but where we generally don't know where to put them. If you look at the [UCSC names](https://github.com/dpryan79/ChromosomeMappings/blob/master/GRCh38_ensembl2UCSC.txt) for these, then you can more easily tell if they have an associated chromosome or not.

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
## genetic variance
- intra-species genetic variance (human) : 0.08% (2M ~ 3M bp)
- inter-species genetic variance (between human and chip) : 0.2%

## DNA lesion
- mutation : a change in the sequence of base pair.
- DNA lesion : an abnormal chemical structure in DNA. a section of a DNA molecule containing a primary damaged site i.e. a base alteration, a base deletion, a sugar alteration or a strand break. Replication before repair, or inefficient repair, can result in the fixation of a primary lesion as a permanent mutation.

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
- Pan-cancer analysis of whole genomes | [paper](https://www.nature.com/articles/s41586-020-1969-6) [docekr](https://dockstore.org/organizations/PCAWG/collections/PCAWG)
  - 91% of tumours had at least one identified driver mutation, with an average of 4.6 drivers per tumour identified, showing extensive ariation across cancer types.
  - For coding point mutations, the average was 2.6 drivers per tumour.
  - Only 13% of driver point mutations were non-coding in PCAWG. Nonetheless, 25% of PCAWG tumours bear at last one putative non-coding driver point mutation, and one third affected the TERT promotoer (9% of PCAWG tumours). Overall, non-coding driver point mutations are less frequent than coding driver mutations.
