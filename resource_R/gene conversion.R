# gene conversion
# ref : https://bioinformatics.stackexchange.com/questions/5229/converting-gene-symbol-to-ensembl-id-in-r

## EnsDb.Hsapiens.v79
BiocManager::install("EnsDb.Hsapiens.v79")
library(EnsDb.Hsapiens.v)

# 1. Convert from ensembl.gene to gene.symbol
ensembl.genes <- c("ENSG00000150676", "ENSG00000099308", "ENSG00000142676", "ENSG00000180776", "ENSG00000108848", "ENSG00000277370", "ENSG00000103811", "ENSG00000101473")
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

# 2. Convert from gene.symbol to ensembl.gene
geneSymbols <-  c('DDX26B','CCDC83',  'MAST3', 'RPL11', 'ZDHHC20',  'LUC7L3',  'SNORD49A',  'CTSH', 'ACOT8')
geneIDs2 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))

## biomart
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
getBM(attributes='hgnc_symbol', 
      filters = 'ensembl_gene_id', 
      values = ensemblsIDS, 
      mart = ensembl)
