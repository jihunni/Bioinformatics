library(readr)
library(tidyverse)
library(pheatmap)
library(ggplot2)

### in vivo ###
#load an example file
GSM465 = read_tsv("../data/inVitro/GSM1169465.out",col_names = FALSE)

#initial parameter for for loop
out = data.frame(geneID=GSM465$X1, symbol = GSM465$X2)
##i = 491
for(i in 465:485){
    read = read_tsv(paste0("../data/inVitro/GSM1169", i, ".out"), col_names = FALSE)
    colnames(read)[3] = paste0("GSM1169", i)
    if(out$gene == read[[1]]){
        out = cbind(out, read[,3])
    }
    else{
        print("ERROR")
    }
}

str(out)
out = out[-c(53701:53705),]
rownames(out) = out$geneID

condition = read_tsv("../data/inVitro/sample_info.txt", col_names = FALSE)
split = str_split_fixed(condition$X2, pattern = "_", n = 3)
condition = cbind(condition, split)
condition = data.frame(sample=colnames(out)[-c(1,2)], condition=c())
colnames(condition) = c("sample", "info", 'type', 'hour', 'rep')


##save(out, condition, file='inVivo_all_2021.05.15.rda')
save(out, condition, file='inVitro_2021.06.06.rda')

### in vitro ###

###heatmap###
#correlation only with protein coding gene
#anntoation biomart version 2 (no digital below 0)
require("biomaRt")
listMarts()
mart <- useMart("ENSEMBL_MART_MOUSE")
listDatasets(mart)
mart <- useDataset("mmc57bl6nj_gene_ensembl", mart)
filter = listFilters(mart)
ensLookup = annotLookup <- out$geneID
annotLookup <- getBM(
    mart=mart,
    attributes=c("arrayexpress", "ensembl_gene_id",
                 "gene_biotype", "external_gene_name"),
    filter="arrayexpress",
    values=ensLookup,
    uniqueRows=TRUE)

annotLookup <- data.frame(
    ensLookup[match(annotLookup$ensembl_gene_id, ensLookup)], annotLookup)

colnames(annotLookup) <- c(
    "original_id",
    c("ensembl_gene_id", "gene_biotype", "external_gene_name"))


out_label = out[,c(1,2)]
out_label = left_join(out_label, annotLookup, by=c("geneID" = "arrayexpress"))
out_anno = cbind(out_label, out[,-c(1,2)])

save(out_anno, file = 'inVitro_all_annotation_2021.06.06.rda')
rm(out_label, split, mart, list, ensLookup, filter, annotLookup, out)
GS_hallamrk = read_tsv('../data/h.all.v7.4.symbols.gmt',col_names = FALSE)
GS_c2.kegg = read_tsv('../data/c2.cp.kegg.v7.4.symbols.gmt',col_names = FALSE)

#normalization
colSums(out_anno[,-c(1:5)])
library(DESeq2)
se = out_anno[,-c(1:5)]
dds <- DESeqDataSetFromMatrix(se, colData = condition, design = ~ 1 )
dds <- estimateSizeFactors( dds )
sizeFactors(dds)
logcounts <- log2( counts(dds, normalized=TRUE) + 1 )
out_anno[,-c(1:5)] = logcounts
rm(logcounts)
#heatmap
#total
out_anno_select = out_anno

#1. select protein coding gene
out_anno_select = out_anno[out_anno$gene_biotype=='protein_coding',]
out_anno_select = drop_na(out_anno_select, gene_biotype)

filter = as.vector(GS_c2.kegg[GS_c2.kegg$X1 == 'KEGG_RIBOSOME',-c(1,2)])
out_anno_select = out_anno[toupper(out_anno$symbol) %in% filter,] #toupper()
out_anno_select = drop_na(out_anno_select, gene_biotype)
#2. select Immunoglobulin (Ig) variable chain and T-cell receptor (TcR) genes imported or annotated according to the IMGT.
# out_anno_select = out_anno[out_anno$gene_biotype !='protein_coding' & out_anno$gene_biotype !='polymorphic_pseudogene' ,]
# out_anno_select = drop_na(out_anno_select, gene_biotype)
# table(is.na(out_anno_select))

#sample correlation
library("RColorBrewer")
correlationMatrix = cor(out_anno_select[,-c(1:5)])
rownames(correlationMatrix) <- condition$type
colnames(correlationMatrix) <- condition$sample
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#annotation
color_code=data.frame(sample = rownames(correlationMatrix))
rownames(color_code) = colnames(correlationMatrix)
#heatmap
pheatmap(correlationMatrix,
         clustering_distance_rows="correlation",
         clustering_distance_cols="correlation",
         # kmeans_k =4,  #the number of kmeans clusters to make
         #  cutree_rows = 2,
         #  cutree_cols = 2,
         col=colors,
         annotation_col = color_code,
         show_colnames = T,
         show_rownames = T,
         display_numbers = T,
         labels_row = condition$info,
         main = "sample correlation (KEGG-Cell cycle)")

#PCA plot - gene 
input_matrix = as.matrix(log10(out_anno_select[,-c(1:5)]+1))
input_matrix = t(input_matrix)
pca.out = prcomp(input_matrix)

pca.result <- data.frame(id = rownames(input_matrix), X1 = pca.out$x[,1], X2 = pca.out$x[,2])
pca.result = left_join(pca.result, condition, by=c("id" = "sample"))
pca.result = pca.result[order(pca.result$id),]
rm(pca.out, input_matrix)

ggplot(pca.result, aes(x = X1, y = X2, color=type, shape=hour)) + 
    geom_point(size = 5) + 
#    scale_color_manual(values=c('green','red', 'blue'))+
   # geom_text(aes(label = id), col = "black") + 
    #geom_abline(intercept = 0, slope = 0) +
    ggtitle("PCA plot for ribosome")
