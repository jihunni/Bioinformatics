# database for gene regulation
## raw data
- ENCODE
- Cistrome DB
- GTRD
## TF
- HOCOMOCO
- TRANSFAC

# Pre-analysis
## Pre-alignment QC
Tool: 
- FastQC
- cutdapt
- AdapterRemoval v2
- Skewer
- trimmomatic


consideration:
- sequence depth
- sequence quality
- GC content
- Duplication
- Length distribution
- K-mer
- Adapter


## Alignment
Tool:
- BEA-MEM
- Bowtie2


considerations:
- Unique mapping %
- duplicated reads %
- fregment length distribution
- GC bias 


## Post alignment processing and QC
Tools:
- SAMtools
- ATACseqQC
- MultiQC

# Core analysis
count-based peak callers:
- MACS2
- HOMER

shaped-based peak callers:

Hidden Markov Model based peak callers:
- HMMRATAC

# Advanced analysis
# Integration with multimoics data
