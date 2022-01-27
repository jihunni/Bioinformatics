package

#library
library(readr)
library(ggplot2)
library(tidyverse)
library(readxl)
library("Hmisc") #correlation matrix with p-value
library(corrplot) #visualizing correlation matrix
library(beepr) #beep() : alarm function
library(DESeq2)
library("apeglm") #shrink log fold change
library(ggraph)
library(igraph)
library("vsn")
library(piano) #report metabolite analysis
#install package
require()

packageVersion("TCGAbiolinks")

install.packages('ggraph')

if (!"BiocManager" %in% rownames(installed.packages()))
    install.packages("dendextend")
BiocManager::install("piano")

sessionInfo()

##Error and trials
ERROR: failed to lock directory ‘C:\R\R-4.0.3\library’ for modifying
Try removing ‘C:\R\R-4.0.3\library/00LOCK
1. install.packages("설치하고 싶은 패키지 이름", dependencies = TRUE, INSTALL_opts = "--no-lock")
2. find the file:‘C:\Program Files\R\R-3.6.1\library/00LOCK’, delete the file"00LOCK" then it works again

###Tidyverse###
# Tidyr
## To split one column into two column by deliminator '/'
new_df = tidyr::separate(data= old_df, col= ColName, sep= '/', into=c("col_name1", "colname2"))
## To unite two column into one column by deliminator '/'
new_df = tidyr::unite(data = old_df, new_colName, old_ColName1, old_ColName2, sep = '/')
# dpyler
##extract into new data.frame
new_data.frame=dpyler::select(data.frame, colum name)
new_data.frame=dpyler::filter(data.frame, boolian for all row)

# data.table ##########################                           
library(data.table)
df = fread('file_name.tsv', sep='\t', stringsAsFactors=TRUE, header=FALSE, col.names=c('name', 'start', 'end', 'peak', 'phyloP', 'gene'))
                              
##statistics
statistics = data.frame %>%
    filter(colname==condition) %>%
    group_by(colnames) %>%
    summarise(number=n()) %>%
    na.omit()
                              
###Bioinformatics###
#convert gene name to ensemble gene ID
library("AnnotationDbi")
library("org.Hs.eg.db")
options(connectionObserver = NULL)
gene_id = mapIds(org.Hs.eg.db,
                 keys=TCGA_LIHC_gene_count$sample, 
                 column="ENSEMBL",
                 keytype="SYMBOL",
                 multiVals="first")
gene_id = data.frame(gene_id)
nrow(gene_id)
length(which(gene_id$gene_id != 'NA')) #17311
TCGA_LIHC_gene_count = cbind(gene_id, TCGA_LIHC_gene_count) #add sample column at first


# to map EntrezID to gene symbol
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#keytypes(TxDb.Hsapiens.UCSC.hg38.knownGene)
#columns(EnsDb.Hsapiens.v86)
entrez = anno$geneId
annotations_edb = AnnotationDbi::select(EnsDb.Hsapiens.v86, keys=entrez, columns=c("GENEID", "GENENAME", "SYMBOL"), keytype="ENTREZID")
annotations_edb$ENTREZID = as.character(annotations_edb$ENTREZID)
  anno = anno %>%
    left_join(annotations_edb, by=c("geneId"="ENTREZID")) 

#annotation, biomart
require("biomaRt")
listMarts()
mart <- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(mart)
mart <- useDataset("hsapiens_gene_ensembl", mart)
filter = listFilters(mart)

ens <- c("ENSG00000100601.5", "ENSG00000178826.6",
         "ENSG00000243663.1", "ENSG00000138231.8")
ensLookup <- gsub("\\.[0-9]*$", "", ens)

annotLookup <- getBM(
    mart=mart,
    attributes=c("ensembl_transcript_id", "ensembl_gene_id",
                 "gene_biotype", "external_gene_name"),
    filter="ensembl_gene_id",
    values=ensLookup,
    uniqueRows=TRUE)

annotLookup <- data.frame(
    ens[match(annotLookup$ensembl_gene_id, ensLookup)],
    annotLookup)

colnames(annotLookup) <- c(
    "original_id",
    c("ensembl_transcript_id", "ensembl_gene_id",
      "gene_biotype", "external_gene_name"))

# BioMart : ensemble_transcript
library(biomaRt)
listMarts()
mart <- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(mart)
mart <- useDataset("hsapiens_gene_ensembl", mart)
filter = listFilters(mart) # to check the type of key

ensLookup = peakAnno$transcriptId
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id", "ensembl_gene_id",
               "gene_biotype", "external_gene_name"),
  filter="ensembl_transcript_id_version",
  values=ensLookup,
  uniqueRows=TRUE)
                              
#anntoation biomart version 2 (no digital below 0)
require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
listMarts()
mart <- useDataset("hsapiens_gene_ensembl", mart)
listDatasets(mart)
filter = listFilters(mart)

ensLookup <- gene_id
annotLookup <- getBM(
  mart=mart,
  attributes=c(,"ensembl_gene_id","gene_biotype", "external_gene_name"),
  filter="ensembl_gene_id",
  values=ensLookup,
  uniqueRows=TRUE)

annotLookup <- data.frame(
  ensLookup[match(annotLookup$ensembl_gene_id, ensLookup)], annotLookup)

colnames(annotLookup) <- c(
  "original_id",
  c("ensembl_transcript_id", "ensembl_gene_id",
    "gene_biotype", "external_gene_name"))
    
# convert Entrez ID into Ensembl
library("biomaRt")
table(gene_dependency_threshold0.001_cytoscape$gene2_entrez  %in% gene_dependency_threshold0.001_cytoscape$gene1_entrez) 
    # one converting list is enough
entrezgene_ = gene_dependency_threshold0.001_cytoscape$gene1_entrez

#listMarts()
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
filter = listFilters(mart)
convert_genes <- getBM(filters="entrezgene_id", 
               attributes=c("ensembl_gene_id","entrezgene_id"), 
               values=entrezgene, mart=mart)

#TCGAutils ######################################
# To convert TCGA UUID to Barcode
library(TCGAutils)
library(curatedTCGAData)
library(MultiAssayExperiment)
library(RTCGAToolbox)
library(BiocFileCache)
library(rtracklayer)
library(R.utils)

UUIDtoBarcode("93d5a265-c116-4179-9fd7-8057d4b53527", from_type = "file_id")

#To convert file UUID into Barcode 
library(tidyverse)
metadata = read_tsv('TCGA-COAD_ATACseq_metadata.tsv')


for (iteration in 1:length(metadata$`File UUID`)){
  file_uuid = metadata$`File UUID`[iteration]
  barcode = UUIDtoBarcode(file_uuid, from_type = "file_id")
  print(barcode)
  metadata$'barcode'[iteration] = barcode$associated_entities.entity_submitter_id
}
save(metadata, file='TCGA-COAD_ATACseq_metadata.rda')

