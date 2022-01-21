library(readr)
library(tidyverse)
library(DESeq2)
library("ggplot2")
library("apeglm") #Library for Log fold change shrinkage
#library("ashr") #Library for Log fold change shrinkage
library("pheatmap")
library("Rtsne")

# data pre-process
load(propensity_normalized) #clinic data with hepatitis, normalized propensity
load(tumor_noHepatitis.propensity.data, file='tumor_noHepatitis.propensity.data_2021.05.06.rda')

#statistics 
tumor_noHepatitis.propensity.data = tumor_noHepatitis.propensity.data[!tumor_noHepatitis.propensity.data$weights == 0,]
table(tumor_noHepatitis.propensity.data$sex)
Female   Male   total
33     32       65

table(substr(colnames(HCC_geneCount_in_Clinic),start=13,stop=15))


HCC_geneCount_in_Clinic = HCC_geneCount[tumor_noHepatitis.propensity.data$case]
colnames(HCC_geneCount_in_Clinic) == tumor_noHepatitis.propensity.data$case
rownames(HCC_geneCount_in_Clinic) = substr(rownames(HCC_geneCount), start=1, stop=15)

clinic_process = tumor_noHepatitis.propensity.data


HCC_geneCount_in_Clinic = 2^HCC_geneCount_in_Clinic - 1 #reverse log2(count+1) transformation

#save input data form DESeq2
save(HCC_geneCount_in_Clinic, clinic_process, file='DESeq2_input_2021.05.07.rda')
rm(HCC_geneCount)



###DEseq###
#to prepare DESeq2 input file  : countData, condition
countData = HCC_geneCount_in_Clinic
countData = data.frame(lapply(countData, as.integer))
rownames(countData) = substr(rownames(HCC_geneCount_in_Clinic),start=1, stop=15)

condition = factor(clinic_process$sex)
weight_vector = clinic_process$weights
rm(tumor_noHepatitis.propensity.data)


###Differential Expression analysis (DESeq2)###
#Did you check rowname of countdata?
#To prepare
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = DataFrame(condition),
                              design = ~ condition) #generate the deseq data set
#filter
dds <- dds[ rowSums(counts(dds)) > 0, ] #remove genes with zero counts

#Note on filter level
dds$condition
#dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
dds$condition <- relevel(dds$condition, ref = "Female") #to specify the reference level
#dds$condition <- droplevels(dds$condition)

#DEA analysis
dds <- estimateSizeFactors(dds)
#dds.normalized = DESeq2::counts(dds, normalized = T) #Accessors for the 'counts' slot of a DESeqDataSet object
#write.table(dds.normalized, "dds.normalized_count.txt", sep="\t")
dds$sizeFactor = dds$sizeFactor * weight_vector
dds = estimateDispersions(dds)
dds = nbinomWaldTest(dds) 

rm(countData, condition, PanCan_HCC_barcode_dropNA, countData,PanCan_HCC_clinic2, HCC_geneCount_dropNA, HCC_geneCount_in_Clinic, weight_vector, clinic_process)

#DEA results
head(sizeFactors(dds))
res = results(dds) #default alpha is 0.1
##summary(res)

#Log fold change shrinkage for visualization and ranking (3 type)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_Male_vs_Female", type="apeglm")


#save the results
save(dds, resLFC, file='dds_result_2021.05.07.rda')


#p-value and adjusted p-value
resOrdered <- resLFC[order(resLFC$pvalue),] #order our results table by the smallest p value
##sum(res$padj < 0.1, na.rm=TRUE) #How many adjusted p-values were less than 0.1?

#if I want to change alpha value into 0.05,
# res05 <- results(dds, alpha=0.05)
# summary(res05)
# sum(res05$padj < 0.05, na.rm=TRUE)


#vsd <- vst(dds, blind = FALSE) #normalization considering tissue 

##MA-plot
plotMA(res, main='MA-plot')
plotMA(res, ylim=c(-7,7), main='MA-plot(male vs female)')
plotMA(resLFC, ylim=c(-10,10), main='MA-plot (male vs female, apeglm)')

#Dispersion Plot
plotDispEsts(dds)

#More information on results columns
mcols(res)$description

#To export results to CSV files
write.csv(as.data.frame(resOrdered), 
          file="PanCan_HCC_male_noHepatitis_vs_female_noHepatitis.csv")
resSig <- subset(resOrdered, padj < 0.1) #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.

#volcano plot
cut_lfc <- 1
cut_pvalue <- 0.01

par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(resLFC)

## Adjusted P values
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot (male tumor without hepatitis vs female without hepatitis, apeglm)", col='grey', cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))

with(subset(topT, padj<cut_pvalue & log2FoldChange>cut_lfc), points(log2FoldChange, -log10(padj), pch=20, col='red', cex=1.5))
with(subset(topT, padj<cut_pvalue & log2FoldChange<(-cut_lfc)), points(log2FoldChange, -log10(padj), pch=20, col='blue', cex=1.5))
## Add lines for FC and P-value cut-off
abline(v=0, col='black', lty=3, lwd=1.0)
abline(v=-cut_lfc, col='black', lty=4, lwd=2.0)
abline(v=cut_lfc, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(topT$padj[topT$padj<cut_pvalue], na.rm=TRUE)), col='black', lty=4, lwd=2.0)

rm(cut_lfc,cut_pvalue,topT, resOrdered)




### Data transformations and visualization ###
# Extracting transformed values
vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE)

# Effects of transformations on the variance
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd)) # to plot row standard deviations versus row means.
head(assay(vsd), 3)

meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

# Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")

rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- colnames(countData)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#annotation
color_code=data.frame(sample = rownames(sampleDistMatrix))
rownames(color_code) = colnames(sampleDistMatrix)
#heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         #  cutree_rows = 2,
         #  cutree_cols = 2,
         col=colors,
         annotation_col = color_code,
         show_colnames = F,
         show_rownames = F,
         main = "Sample-to-sample distances (HCC, no hepatitis)")

# Principal component plot of the samples
##plotPCA(vsd, intgroup="condition")

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    ggtitle("PCA plot_PanCan(HCC, no hepatitis)")+
    coord_fixed()

#tSNE
set.seed(1) # for reproducibility
tsne <- Rtsne(sampleDistMatrix, dims = 2, perplexity=20
              , verbose=TRUE, max_iter = 500)

##tSNE visualizing
colors = rainbow(dim(sampleDistMatrix)[1])
names(colors) = dds$condition
par(mgp=c(2.5,1,0))

color_annotation = gsub("Female","red", dds$condition)
color_annotation = gsub("Male","blue", color_annotation)
)

plot(tsne$Y, col=color_annotation, main="tSNE_PanCan(HCC, no hepatitis)", xlab = "tSNE1", ylab = "tSNE2", pch=16 )
legend("topleft", 
       c("Male", "Female"),
       col=c("blue", "red"),
       pch=c(16, 16))

unique(condition)

rm(sampleDistMatrix, sampleDists, percentVar, colors, color_code, )
rm(sampleDists, colors, percentVar, sampleDistMatrix, color_code, pcaData,color_annotation)
rm(tsne, vsd, ntd)

#extract resLFC list
library(tidyverse)
df = data.frame(resLFC)
df = drop_na(df)
df_up = df[df$log2FoldChange > 1 & df$padj < 0.05,]
df_down = df[df$log2FoldChange < -1 & df$padj < 0.05,]
write.csv(df, 
          file="hepatitis_total_2021.05.07.csv")
write.csv(df_up, 
          file="hepatitis_up_2021.05.07.csv")
write.csv(df_down, 
          file="hepatitis_down_2021.05.07.csv")

save(df, df_up, df_down, file='dataframe_resLFC_dropNA_2021.05.07.rda')
rm(dds, res, resLFC, male)

#Prepare input file for DAVID
##foreground gene
load(dataframe_resLFC_dropNA)
foreground_geneList = substr(rownames(male), start=1, stop=15)
write(foreground_geneList, 'David_propensity_foreground_geneList_2021.04.27_male.txt')

##background gene
load(dataframe_resLFC)
background_geneList = substr(rownames(dataframe_resLFC), start=1, stop=15)
background_geneList = background_geneList[!background_geneList %in% foreground_geneList] 
write(background_geneList, 'David_propensity_background_geneList_2021.04.27_male.txt')

#Input file for GSEA
rm(description,geneID_ensemble )
input_GSEA_normalized_count = counts(dds, normalized =TRUE)
rownames(input_GSEA_normalized_count) = substr(rownames(input_GSEA_normalized_count), start=1, stop=15)
description = data.frame(DESCRIPTION = rep(NA, nrow(input_GSEA_normalized_count)))
geneID_ensemble = data.frame(NAME = rownames(input_GSEA_normalized_count))
input_GSEA_normalized_count = cbind(geneID_ensemble, description, input_GSEA_normalized_count)
write.table(input_GSEA_normalized_count, "GSEA input file_normalizedCount_ensembleID_2021.05.07.txt",quote=F,sep="\t", row.names = F)

rm(description,geneID_ensemble )



##prepare the cls file for GSEA
library(plyr)
#cls = revalue(dds$condition, c("Female" = 0, "Male" = 1)) #change name of factor
cls = as.character(dds$condition)

unique(dds$condition)
table(dds$condition)
unique(cls)
cat("65 2 1", "\n", #num_sample ; num_class ; always 1
    "# Female Male", "\n", #The order of the labels on the third line determines the association of class names and class labels
    cls, file = "GSEA_input_2021.05.07.cls", append = FALSE)

rm(cls, input_GSEA_normalized_count)


#You need to remove space and check the seqeunce in GSEA program
rm(background_geneList, foreground_geneList, cls, dataframe_resLFC, dataframe_resLFC_dropNA)






#compare tumor and tumor Free
tumor_male = read_csv('./output_tumor_2021.04.26/PanCan_HCC_male_2021.04.26.csv')
tumor_male$gene = substr(tumor_male$gene,start=1,stop=15)
tumorFree_male = read_csv('./output_tumorFree_2021.04.26/PanCan_HCC_male_2021.04.26.csv')

table(tumor_male$gene %in% tumorFree_male$gene)
tumor_male_only = tumor_male[!tumor_male$gene %in% tumorFree_male$gene, ]



#convert gene name to ensemble gene ID
options(connectionObserver = NULL)
library("AnnotationDbi")
library("org.Hs.eg.db")

keytypes(org.Hs.eg.db)
gene_id = mapIds(org.Hs.eg.db,
                 keys=rownames(df_B_up), 
                 column="ENSEMBL",
                 keytype="SYMBOL",
                 multiVals="first")
gene_id = data.frame(gene_id)
nrow(gene_id)
length(which(gene_id$gene_id != 'NA')) #17311
TCGA_LIHC_gene_count = cbind(gene_id, TCGA_LIHC_gene_count)

#annotation
require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

#ens <- c("ENSG00000100601.5", "ENSG00000178826.6",
#         "ENSG00000243663.1", "ENSG00000138231.8")
#ensLookup <- gsub("\\.[0-9]*$", "", ens)

ens = ensLookup = rownames(df_down)

annotLookup <- getBM(
    mart=mart,
    attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name"),
    filter="ensembl_gene_id",
    values=ensLookup,
    uniqueRows=TRUE)

annotLookup <- data.frame(
    ens[match(annotLookup$ensembl_gene_id, ensLookup)],
    annotLookup)

colnames(annotLookup) <- c(
    "original_id",
    c("ensembl_gene_id", "gene_biotype", "external_gene_name"))

##write.csv(as.data.frame(annotLookup), 
          file="hepatitis_total_annotatio.csv")
write.table(as.data.frame(annotLookup), file = "hepatitis_down_annotation.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote=FALSE, append=FALSE)


rm(mart, ens, ensLookup, annotLookup)

save(tumor_male, tumorFree_male, tumor_male_only, file='tumor_and_tumorFree_list.rda')
rm(tumor_male, tumorFree_male, tumor_male_only)
rm(ens, ensLookup, mart)

#visualize pathway
library(pathview)
##data(gse16873.d)
##data(demo.paths)
##data(gene.idtype.list)
##demo.paths$sel.paths
pathview_data_lfc = df$log2FoldChange
pathview_data_padj = df$padj
names(pathview_data_lfc)  = names(pathview_data_padj) = rownames(df)
input_pathway.id = "04310"
input_out.suffix = "Wnt signaling pathway"
pathview(gene.data = pathview_data_lfc, gene.idtype ='ENSEMBL', pathway.id = input_pathway.id, species = "human", kegg.dir="C:/R/default_working_directory/Data/KEGG_pathway/",  limit=list(gene=1, cpd=10), out.suffix = paste0(input_out.suffix,'_lfc_low'))
pathview(gene.data = pathview_data_lfc, gene.idtype ='ENSEMBL', pathway.id =input_pathway.id, species = "human", kegg.dir="C:/R/default_working_directory/Data/KEGG_pathway/",  limit=list(gene=5, cpd=10), out.suffix = paste0(input_out.suffix,'_lfc_high'))

pathview(gene.data = pathview_data_padj, gene.idtype ='ENSEMBL', pathway.id = input_pathway.id, species = "human", kegg.dir="C:/R/default_working_directory/Data/KEGG_pathway/",  limit=list(gene=1, cpd=10), out.suffix = paste0(input_out.suffix,'_adjp_low'))
pathview(gene.data = pathview_data_padj, gene.idtype ='ENSEMBL', pathway.id =input_pathway.id, species = "human", kegg.dir="C:/R/default_working_directory/Data/KEGG_pathway/",  limit=list(gene=5, cpd=10), out.suffix = paste0(input_out.suffix,'_adjp_high'))
rm(input_pathway.id, input_out.suffix)

####clustering - 2021.05.19###
library(tidyverse)
library(pheatmap)
library(ggplot2)

normalized_Genecount = counts(dds, normalized =TRUE)
colnames(normalized_Genecount) = gsub("[.]","-", colnames(normalized_Genecount))
save(normalized_Genecount, file='normalized_geneCount.rda')

#read annotation
output_2d = read_tsv('./2d_differential_gene_GSanno_2021.05.17.txt')
output_2d = unique(output_2d) #remove duplication!!
#load data
condition = data.frame(sample=colnames(normalized_Genecount), sex = clinic_process$sex)
condition = condition[order(condition$sample),]


#filter
filter_all = data.frame(gene = rownames(normalized_Genecount))
filter_all = left_join(filter_all, output_2d, by=c("gene" = "gene")) #data.frame
##filter_gene = filter_all$gene_biotype == 'lncRNA' #boolian, 
table = table(filter_all$GS_c2.kegg)
##filter_gene = filter_all$GS_hallmark == 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY'
filter_gene = filter_all$GS_c2.kegg == 'KEGG_GLUTATHIONE_METABOLISM' #boolian, 
##filter_gene = filter_all$gene_biotype == 'lncRNA' & (filter_all$padj_sex <= 0.05 & filter_all$padj_normal <= 0.05) & (filter_all$LFC_sex > 1 | filter_all$LFC_sex < -1) & (filter_all$LFC_normal > 1 | filter_all$LFC_normal < -1) 
##filter_gene = (filter_all$padj_sex <= 0.05 & filter_all$padj_normal <= 0.05) & (filter_all$LFC_sex > 1 | filter_all$LFC_sex < -1) & (filter_all$LFC_normal > 1 | filter_all$LFC_normal < -1) 
filter_gene[is.na(filter_gene)] = FALSE
table(filter_gene)

#correlation matrix
library(pheatmap)
library("RColorBrewer")
correlationMatrix = cor(normalized_Genecount[filter_gene,])
rownames(correlationMatrix) <- condition$sex
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
         display_numbers = F,
         labels_row = condition$info,
         main = "HCC correlation (sex-specific coding gene)")

ggsave(paste0("./figure/KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE",".png"), width=10, height = 7, units='cm', limitsize = FALSE)

#PCA plot
##pca.out = prcomp(t(as.matrix(log10(normalized_Genecount[filter_gene,]+1), cen=T,sca=T))) #gene PCA plot
pca.out = prcomp(t(as.matrix(log10(normalized_Genecount[filter_gene,]+1), cen=T,sca=T))) #sample PCA plot
    ##log(x+1) transformation to remove x=0 case (inf.)
color_annotation = gsub("Female","red", dds$condition)
color_annotation = gsub("Male","blue", color_annotation)

plot(pca.out$x[,1:2], col=color_annotation, pch=16, main='PCA plot (sex-specific gene)')
legend("topleft", # "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center"
       c("Male", "Female"),
       col=c("blue", "red"),
       pch=c(16, 16))

pca.out.df = data.frame(sample=colnames(normalized_Genecount))
pca.out.df2 = data.frame(sample=colnames(normalized_Genecount), PC1 = pca.out$x[,1], PC2= pca.out$x[,2])
pca.out.df = left_join(pca.out.df, pca.out.df2, by=c("sample" = "sample"))
pca.out.df = left_join(pca.out.df, condition, by=c("sample" = "sample"))
rm(pca.out.df2)

ggplot(pca.out.df, aes(PC1, PC2, color=sex)) +
    geom_point(size=3) +
    xlab(paste0("PC1")) +
    ylab(paste0("PC2")) +
    ggtitle("PCA plot - HCC (sex-specific lncRNA)")+
    coord_fixed()

save(filter_all, normalized_Genecount, condition, file='input_clusting_2021.05.19.rda')

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
