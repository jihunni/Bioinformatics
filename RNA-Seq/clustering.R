####clustering - 2021.05.19###
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggplot2)

###normalization###
#data preprocessing for log count normalization
HBV = HCC_clinic4$HBV
HCV = HCC_clinic4$HCV
condition = rep(NA, nrow(HCC_clinic4))
names(condition) = HCC_clinic4$sample

for (i in 1:nrow(HCC_clinic4)){
    if(HBV[i] == FALSE & HCV[i] == FALSE){
        condition[i] = 'no'
    }
    if(HBV[i] == TRUE & HCV[i] == FALSE){
        condition[i] = 'HBV'
    }
    if(HBV[i] == FALSE & HCV[i] == TRUE){
        condition[i] = 'HCV'
    }    
        
    if(HBV[i] == TRUE & HCV[i] == TRUE){
        condition[i] = 'both'
    }    
}
check = cbind(HCC_clinic4, condition)
rm(check, HBV, HCV, i)
condition = cbind(sample = rownames(condition), condition)
condition = data.frame(infection=condition)
condition = cbind(condition, sex = HCC_clinic4$gender)

condition['tumor'] = substr(rownames(condition), start=14, stop=15)
condition$tumor = gsub("01", "tumor", condition$tumor)
condition$tumor = gsub("02", "tumor", condition$tumor)
condition$tumor = gsub("11", "normal", condition$tumor)

#log count normalization
#HCC_geneCount_in_Clinic = 2^HCC_geneCount_in_Clinic - 1
se = HCC_geneCount #data type : list
se = 2^se - 1
se = data.frame(lapply(se, as.integer))

dds <- DESeqDataSetFromMatrix(se, colData = condition, design = ~ 1 )
dds <- estimateSizeFactors( dds )
sizeFactors(dds)
logcounts <- log2(counts(dds, normalized=TRUE) + 1)
rownames(logcounts) = substr(rownames(HCC_geneCount),start=1,stop=15)
rm(se, dds, HCC_geneCount)

###if we have already normlaized datafrom DESeq2
# normalized_Genecount = counts(dds, normalized =TRUE)
# colnames(normalized_Genecount) = gsub("[.]","-", colnames(normalized_Genecount))

# #read annotation
# output_2d = read_tsv('./2d_differential_gene_GSanno_2021.05.17.txt')
# output_2d = unique(output_2d) #remove duplication!!
# #load data
# condition = data.frame(sample=colnames(normalized_Genecount), sex = clinic_process$sex)
# condition = condition[order(condition$sample),]

#2nd data set :: noH with propensity
noH_normalized = counts(dds, normalized =TRUE)
colnames(noH_normalized) = gsub("[.]","-",colnames(noH_normalized))
noH_normalized <- log2(noH_normalized + 1)


#annotation
require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
listMarts()
mart <- useDataset("hsapiens_gene_ensembl", mart)
listDatasets(mart)
filter = listFilters(mart)

ensLookup <- rownames(noH_normalized)
annotLookup <- getBM(
    mart=mart,
    attributes=c("ensembl_gene_id","gene_biotype", "external_gene_name"),
    filter="ensembl_gene_id",
    values=ensLookup,
    uniqueRows=TRUE)

##gene set inclusion
load('../Gene Set Enrichment Analysis/Gene Set inclusion/GS_merge_hallmark_2021.05.17.rda')
GS = drop_na(GS)
annotLookup = left_join(annotLookup, GS, by=c("external_gene_name"="gene"))
colnames(annotLookup)[4] = 'hallmark'

load('../Gene Set Enrichment Analysis/Gene Set inclusion/GS_merge_c2.kegg_2021.05.17.rda')
GS = drop_na(GS)
annotLookup = left_join(annotLookup, GS, by=c("external_gene_name"="gene"))
colnames(annotLookup)[5] = 'kegg'


##merge with log count
library(tidyverse)
normalized_Genecount = as.data.frame(logcounts)
normalized_Genecount = cbind(geneID = rownames(normalized_Genecount), normalized_Genecount)
normalized_Genecount = left_join(annotLookup, normalized_Genecount, by=c("ensembl_gene_id"="geneID"))
rownames(normalized_Genecount) = normalized_Genecount$ensembl_gene_id
save(normalized_Genecount, condition, file='normalized_geneCount_2021.06.12.rda')

#merge with noH_normalized
library(tidyverse) 
noH_normalized = as.data.frame(noH_normalized)
noH_normalized = cbind(geneID = rownames(noH_normalized), noH_normalized)
noH_normalized = left_join(annotLookup, noH_normalized, by=c("ensembl_gene_id"="geneID"))
rownames(noH_normalized) = noH_normalized$ensembl_gene_id
save(noH_normalized, condition, file='noH_normalized_2021.06.12.rda')

#filter
normalized_Genecount = noH_normalized
##colnames(normalized_Genecount)[-c(1:5)] == condition$sample
filter_sample = condition$tumor == 'tumor' & condition$infection == 'no'

filter_gene = normalized_Genecount$gene_biotype == 'protein_coding' #boolian, 
##filter_gene = normalized_Genecount$GS_hallmark == 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY'
filter_gene = normalized_Genecount$kegg == 'KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS'
filter_gene = normalized_Genecount$kegg == 'KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS' | normalized_Genecount$kegg == 'KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450'| normalized_Genecount$kegg == 'KEGG_STEROID_HORMONE_BIOSYNTHESIS'|normalized_Genecount$kegg == 'KEGG_RIBOSOME'|normalized_Genecount$kegg == 'KEGG_PORPHYRIN_AND_CHLOROPHYLL_METABOLISM'|normalized_Genecount$kegg == '	KEGG_RETINOL_METABOLISM'|normalized_Genecount$kegg == 'KEGG_PRIMARY_BILE_ACID_BIOSYNTHESIS'|normalized_Genecount$kegg == 'KEGG_DRUG_METABOLISM_CYTOCHROME_P450'

##filter_gene = normalized_Genecount$gene_biotype == 'lncRNA' & (normalized_Genecount$padj_sex <= 0.05 & normalized_Genecount$padj_normal <= 0.05) & (normalized_Genecount$LFC_sex > 1 | normalized_Genecount$LFC_sex < -1) & (normalized_Genecount$LFC_normal > 1 | normalized_Genecount$LFC_normal < -1) 
##filter_gene = (normalized_Genecount$padj_sex <= 0.05 & normalized_Genecount$padj_normal <= 0.05) & (normalized_Genecount$LFC_sex > 1 | normalized_Genecount$LFC_sex < -1) & (normalized_Genecount$LFC_normal > 1 | normalized_Genecount$LFC_normal < -1) 
filter_gene[is.na(filter_gene)] = FALSE
table(filter_gene)

input_matrix = normalized_Genecount[filter_gene, c(rep(FALSE, 5), filter_sample)]
rownames(input_matrix) = rownames(normalized_Genecount)[filter_gene]

input_matrix = noH_normalized[filter_gene,-c(1:5)]

#PCA plot
##library(ggplot2)
##pca.out = prcomp(t(as.matrix(log10(normalized_Genecount[filter_gene,]+1), cen=T,sca=T))) #gene PCA plot
pca.out = prcomp(t(as.matrix(input_matrix[rowSums(input_matrix)>0,])), cen=T,sca=T)
pca.result <- data.frame(id = colnames(input_matrix), X1 = pca.out$x[,1], X2 = pca.out$x[,2])
pca.result = left_join(pca.result, condition, by=c("id" = "sample"))
pca.result = pca.result[order(pca.result$id),]
rm(pca.out)

ggplot(pca.result, aes(X1, X2, color=sex, shape = infection)) +
    geom_point(size=2) +
    xlab(paste0("PC1")) +
    ylab(paste0("PC2")) +
    ggtitle("PCA plot")

ggplot(pca.result, aes(X1, X2, color=infection, shape = sex)) +
    geom_point(size=2) +
    xlab(paste0("PC1")) +
    ylab(paste0("PC2")) +
    ggtitle("PCA plot")

#tSNE
#library(Rtsne)
#set.seed(1) # for reproducibility
tsne <- Rtsne(t(as.matrix(input_matrix[rowSums(input_matrix)>0,])), dims = 2, perplexity=30
              , verbose=TRUE, max_iter = 500)

tsne.result <- data.frame(id = colnames(input_matrix), X1 = tsne$Y[,1], X2 = tsne$Y[,2])
tsne.result = left_join(tsne.result, condition, by=c("id" = "sample"))
tsne.result = tsne.result[order(tsne.result$id),]
rm(tsne)

ggplot(tsne.result, aes(X1, X2, color=sex, shape = infection)) +
    geom_point(size=2) +
    xlab(paste0("PC1")) +
    ylab(paste0("PC2")) +
    ggtitle("tSNE plot (male-specific gene set)")

ggplot(tsne.result, aes(X1, X2, color=infection, shape = sex)) +
    geom_point(size=2) +
    xlab(paste0("PC1")) +
    ylab(paste0("PC2")) +
    ggtitle("tSNE plot")

#correlation matrix
library(pheatmap)
library("RColorBrewer")
correlationMatrix = cor(input_matrix)
colnames(correlationMatrix) = rownames(correlationMatrix) <- colnames(input_matrix)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#annotation
color_code=condition

#heatmap
pheatmap(correlationMatrix,
         clustering_distance_rows="correlation",
         clustering_distance_cols="correlation",
         # kmeans_k =4,  #the number of kmeans clusters to make
         #  cutree_rows = 2,
         #  cutree_cols = 2,
         col=colors,
         annotation_col = color_code,
         show_colnames = F,
         show_rownames = F,
         display_numbers = F,
         #labels_row = condition$info,
         main = "HCC correlation (sex-specific coding gene)")

ggsave(paste0("./figure/KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE",".png"), width=10, height = 7, units='cm', limitsize = FALSE)


###hierarcy clustering###
library(dendextend)
distance = as.dist(1 - correlationMatrix)
hclusters = hclust(distance,  method="complete") #“complete linkage” method
dend = as.dendrogram(hclusters)

#color annotation
color_hclusters = as.vector(condition$sex)
names(color_hclusters) = as.vector(condition$sample)
color_hclusters = gsub("Male","blue", color_hclusters)
color_hclusters = gsub("Female","red", color_hclusters)
dend = set(dend, "labels_col", color_hclusters)

plot(dend)


hclusters = color_labels(as.dendrogram(hclusters), labels = names(condition), col=color_hclusters)
plot(hclusters, col = color_hclusters)

sub_grp <- cutree(hclusters, k = 2)
