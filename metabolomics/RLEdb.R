# load human RLEdb
library(readr)
library(tidyverse)

RLEdb_human = read_tsv('RLEs_Human.txt', col_names = FALSE)
colnames(RLEdb_human) = c("gene_symbol", "sub_symbol")

# BioMart : to convert gene symbol to ensembl
library(biomaRt)
listMarts()
mart <- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(mart)
mart <- useDataset("hsapiens_gene_ensembl", mart)
filter = listFilters(mart) # to check the type of key

ensLookup = RLEdb_human$gene_symbol
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id","entrezgene_id","gene_biotype", "external_gene_name"),
  filter="hgnc_symbol",
  values=ensLookup,
  uniqueRows=TRUE)

# HMR subsystem
load('~/R/GEM/genesByHMRSubsystems_df.RData')
RLEdb_human_df = left_join(x= annotLookup, y= genesByHMRSubsystems_df, by=c("ensembl_gene_id"="gene"))
table(RLEdb_human_df$subsystem)

rm(annotLookup, RLEdb_human, filter, mart, genesByHMRSubsystems_df, ensLookup)
save(RLEdb_human_df, file='RLEdb_human_df.Rdata')
write.table(RLEdb_human_df, file = "RLEdb_human_df.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote=FALSE, append=FALSE)
