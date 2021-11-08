# Tumor microenvironment
## paper
- Spatially confined sub-tumor mictroenvironments in pancreatic cancer | Cell, 2021
	- pancreatic ductal adenocarcinoma (PDAC) regional heterogeneity stems from sub-tumor microenvironment (subTMEs)
		- deserted region : thin, spindle-shaped fibroblasts, loose matured fibers, and often keloid or myxoid features
		- reactive region : plump fibroblasts with enlarged nuclei, few acellular components, often rich in inflammatory infiltration
		- intermediate regon
		- sub-tumor microenvironment : co-occurence of deserted region and reactive regions within the same tumor in a spatially confined maner.
	- subTMEs exhibit distinct immune phenotypes and CAF differentiation states
		- deserted region : grow faster
		- reactive region : more motile
	- SubTMEs execute distinct tumor-promoting (reactive) and chemoprotective (deserted) functions
		- deserted region : chemoprotective (chemoresistance)
		- reactive region : proliferative, basal-like, poorly differentiated tumor cell phenotype
	- Intratumoral subTME co-occurence links stromal heterogeneity to patient outcome.
	- method : laser catpure microdissection (LCM), RNA sequencing, shotgu proteomics, single-cell RNA-seq, patient-derived organoids
- 
# Cancer evolution
## medical terminology
- Prevalence is the proportion of a population who have a specific characteristic in a given time period. (e.g. the prevalence of somatic mutations in each sample)
- 
## terminology : genome instability and cancer evolution
- genome stability : 
  - microsatellite instability (MSI) : mutations/epigenetic silencing of mismatch repair gene -> DNA mismatch repair X -> somatic single nucleotide variants (SNVs), small indels
  - chromosomal instability (CIN) : aneuploid, accumulate somatic strctural variants (SVs ; genomic rearrangements) -> copy number variations (CNVs)
  - **Copy neutral loss of heterozygosity** describes a phenomenon whereby one of two homologous chromosomal regions is lost, but various mechanisms have ensured the presence of two identical copies of such region in the genome. As a result, the karyotype appears normal or ‘copy neutral.’
- mutation:
  - synonymous mutation : Synonymous mutations are point mutations that only changes one base pair in the RNA copy of the DNA. They have no real role in the evolution of species since the gene or protein is not changed in any way. Synonymous mutations are actually fairly common, but since they have no effect, then they are not noticed.
  - nonsynonymous mutation : Synonymous mutations are point mutations that changes both one base pair in the RNA copy of the DNA and protein sequence.
- hypodiploid
- hyperdiploid
- clonal
- molecular time
## terminology : sequencing
- Base call quality score (Q score) : Q score is a PHRED-scaled probability ranging from 0-20 inversely proportional to the probability that a single sequenced base is correct. For example, a T base call with Q of 20 is considered likely correct with a confidence P-value of 0.01. Any base call with Q<20 should be considered low quality, and any variant identified where a substantial proportion of reads supporting the variant are of low quality should be considered potentially false positive.
- Read depth : Read depth (or coverage, conventionally a number followed by "×") is the number of independent reads with overlapping
alignment at a locus of interest.
- Variant read number (variant reads): Variant reads is the number of independent sequence reads supporting the presence of a variant. Due to the high error rate of NGS at the per-base call level, calls supported by fewer than 5 variant reads are typically considered to be likely false positive calls.
- variant allele frequency (VAF) : the percentage of sequence reads observed matching a specific DNA variant divided by the overall coverage at that locus. 
  - NGS provides a near random samples
  - heterozygous loci sholud be near 50%
  - homozygous loci sholud be near 100%
  - referene loci sholud be near 0%
  - deviations from these three expected values sholud be considered suspicious as potential errors due to incorrect base calls or alignment.
- Variant quality scores (QUAL) : transformed log-scaled (PHRED) values where, for example, a score of 90 supports the variant call with a P-value of 1×10 -9.
- allele
  - major allele: the most common allele for a given SNP
  - minor allele: the less common allele for a SNP. The MAF is therefore the minor allele frequence. This measure can be used to get a rough idea of the variation of genotypes for a given SNP in a given population, in other words it tells you how common this SNP is.
  - risk allele: in the context of a disease, this is the allele that confers a risk of developing the disease. Most of the time, risk allele = minor allele, as most people will not carry the risk allele. However, in some case, the risk allele can in fact be the major allele.
  - effect allele
  - reference allele
  - wildtype allele
- mutant allele frequency (MAF)
- minor allel (B-allel) frequency (IBAF)
- Tumor purity: the percentage of cancer cells in a solid tumor sample
- copy ratio plot
- 

## tool
- BWA
- Picard, GATK
- MuTect, VarScan2, MuSE : somatic SNV
- VarScan2, GATK, Baylor : Indel
- Sequenza : sample purity, ploidy
## paper
- Heterozygous mutations cause genetic instability in a yeast model of cancer evolution
- A river model to map convergent cancer evolution and guide therapy in RCC
- A renewed model of pancreatic cancer evolution based on genomic rearrangement patterns
- Modeling colorectal cancer evolution
- The clonal evolution of metastatic colorectal cancer
- Genomic evolution and diverse models of systemic metastases in colorectal cancer
- Integrated Multiregional Analysis Proposing a New Model of Colorectal Cancer Evolution 
- `review` Clinical implications of intratumor heterogeneity: challenges and opportunities
- `review` Clonal evolution in cancer | Nature, 2012
- `review` Clonal Heterogeneity and Tumor Evolution: Past, Present, and the Future | Cell, 2017
- `reivew` Integrating genetic and non-genetic determinants of cancer evolution by single-cell multi-omics | Nature
- Temporal dissection of tumorigenesis in primary cancers | Cancer Discovery, 2011 | [paper](https://cancerdiscovery.aacrjournals.org/content/candisc/1/2/137.full.pdf)
- The Life History of 21 Breast Cancers | Cell, 2012
- `method` Estimation of rearrangement phylogeny for cancer genomes | Genome Res. 2012
- Pan-cancer analysis of whole genomes | Nature, 2020 | [paper](https://www.nature.com/articles/s41586-020-1969-6#Fig5)
  - On average, cancer genomes contained 4~5 driver mutations when combining coding and non-coding genomic elements. However, in around 5% of cases no drivers were identified.
  - Chromothripsis, in which many clustered structural variants arise in a single catastrophic event, is frequently an early event in tumour evolution.
- The evolutionary history of 2,658 cancers | Nature, 2020 | [paper](https://www.nature.com/articles/s41586-019-1907-7)
- `review` Modeling colorectal cancer evolution
- `review` Current practices and guidelines for clinical next-generation sequencing oncology testing | Cancer Biol Med 2016 | [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4850126/pdf/cbm-13-1-3.pdf)
  - basic terminology for RNA-seq and mutation  
