# papers
- The chromatin accessibility landscape of primary human cancers | science


# CAGE (cap analysis of gene expression)
to identify and monitor the activity (transcription initiation frequency) of **transcription start sites (TSS)** at single base-pair resultion acorss the genome

# MNase-seq (micrococal nuclease digestion followed by sequencing)
to measure **nucleosome occupancy** in the human genome  
![image](https://user-images.githubusercontent.com/48517782/133800009-dcabf1ed-e5e2-46b5-9891-c6f41f7bdae7.png)

# DNase
- Global mapping of protein-DNA interactions in vivo by digital genomic footprinting | Nature Methods, 2009 | [paper](https://www.nature.com/articles/nmeth.1313)
- An expansive human regulatory lexicon encoded in transcription factor footprints | Nature 2012 | [paper](https://www.nature.com/articles/nature11212) 
-  Coupling transcription factor occupancy to nucleosome architecture with DNase-FLAS
  - refinement of Dna
- Refined DNase-seq protocol and data analysis reveals intrinsic bias in transcription factor footprint identification | Nature Methods, 2013| [paper](https://www.nature.com/articles/nmeth.2762)

# ATAC-seq
![image](https://user-images.githubusercontent.com/48517782/133800455-17b98397-d92d-4f9e-9f7b-d890c6d23c59.png)  

- `method` Rapid, low-input, low-bias construction of shotgun fragment libraries by high-density in vitro transposition | Genome Biology, 2010 | [paper](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-12-r119)
  - Tn5 transpose
  - tagmentation
- `method` Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position | 2013 Nature Methods | [link](https://www.nature.com/articles/nmeth.2688.pdf)
  - chromatin
    - closed : inactive
    - open : active (bounded by TF --> transcription)
  - ATAC-seq
    1. Tn5 transposase gain assess to open/accessible region of genome
    2. tagmentation : transposase insert known DNA seqences (adaptor) into open region
    3. libraries are amplified and sequenced
  - insert size disclose nucleosome position
  - footprint (a deep notch of ATAC-seq) : the presence of DNA-binding protein at each site
- `news` The genome shows its sensitive side | Nature Methods, 2014 | 
- `method` An improved ATACseq protocol reduces background and enables interrogation of frozen tissues | 2017, Nature Method | [link](https://www.nature.com/articles/nmeth.4396.pdf)
- `review` Identifying and mitigating bias in next-generation sequencing methods for chromatin biology | 2014 Nature Reviews Genetics | [paper](https://www.nature.com/articles/nrg3788)
- `review` From reads to insight: a hitchhiker’s guide to ATAC-seq data analysis | 2020, Genome Biology | [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1929-3)
- file format | [link](http://genome.ucsc.edu/FAQ/FAQformat.html#format1)
