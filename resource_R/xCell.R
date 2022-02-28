library(tidyverse)
library(ggplot2)
library("ggpubr") #ggplot with p-value

#load files
load(HCC_RSEM_gene_TPM)
load(tumor_noHepatitis.propensity.data)
clinic_data = rbind(tumor_noHepatitis.propensity.data, tumorFree_noHepatitis.propensity.data)

rownames = rownames(HCC_tpm)
HCC_tpm = HCC_tpm[clinic_data$case]
rownames(HCC_tpm) = rownames
rownames(HCC_tpm) = substr(rownames(HCC_tpm), start=1, stop=15)
save(HCC_tpm, file='tumor_noHepatitis_tpm2_2021.05.07.rda')

#annotation, biomart
require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
ens <- rownames(HCC_tpm)
ensLookup <- gsub("\\.[0-9]*$", "", ens)
ensLookup

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

annotLookup

colnames(annotLookup)[1] ='gene_id'
gene_list = left_join(gene_list,annotLookup, by=c('X1'='gene_id'))
rm(ens, ensLookup,mart)


#merge with count file and annotation file
table(is.na(annotLookup$external_gene_name))
output = HCC_tpm[annotLookup$ensembl_gene_id,]
table(duplicated(annotLookup$external_gene_name))
output_rownames = annotLookup$external_gene_name[!is.na(annotLookup$external_gene_name)]
output = cbind(output_rownames, output)

write.table(output, file = "xCell_input2_PanCan_HCC_2021.05.07.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote=FALSE, append=FALSE) 

#propensity algorithm
library(dplyr)
colnames(output)[-c(1)] == clinic_data$case
weight_vector = data.frame(case=colnames(output)[-1])
vector = tumor_noHepatitis.propensity.data[,c(1,12)]
weight_vector = left_join(weight_vector, vector, by=c("case"="case"))
weight_vector = as.vector(weight_vector)

output_propensity = data.frame(mapply('*',output[,-1],weight_vector$weights))
output_propensity = cbind(output_rownames, output_propensity)


write.table(output_propensity, file = "xCell_input2_PanCan_HCC_propensity_2021.05.07.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote=FALSE, append=FALSE) 

save(output, output_propensity, tumor_noHepatitis.propensity.data, file='2_2021.05.07.rda')
rm(tcga_tpm, vector, weight_vector, example_input, HCC_tpm)
rm(annotLookup, output_rownames)

#read the output from xCell
library(readr)
xCell_result = read_tsv('xCell_output2_PanCan_HCC_propensity_2021.05.07_results.txt')

colnames(xCell_result) = gsub("[.]","-", colnames(xCell_result)) #change colname
colnames(xCell_result)[1] = 'cell_type'
colnames(xCell_pval) = colnames(xCell_result)
xCell_result= as.numeric(xCell_result[,-1])
ex = xCell_result[,-1]
ex = as.numeric(ex)
save(xCell_result, xCell_pval, file='xcell_output_2021.05.07.rda') #save output file


### plot ###
detach("tidyverse")
#grouping in propensity normalized file
#Tumor and TumorFree in Male
male_tumor = clinic_data[clinic_data$sex == 'Male' & clinic_data$status == 'With Tumor',]$case
male_tumorFree = clinic_data[clinic_data$sex == 'Male' & clinic_data$status == 'Tumor Free',]$case
female_tumor = clinic_data[clinic_data$sex == 'Female' & clinic_data$status == 'With Tumor',]$case
female_tumorFree = clinic_data[clinic_data$sex == 'Female' & clinic_data$status == 'Tumor Free',]$case

# Violin plots with box plots inside - for loop
##detach("package:dplyr", unload=TRUE)
for(i in xCell_result$cell_type){
#  i = 'Basophils'
    print(i)
    gene = i
    index = which(xCell_result$cell_type == gene)

    gene_male_tumor = as.numeric(xCell_result[index, colnames(xCell_result) %in% male_tumor])
    gene_male_tumorFree = xCell_result[, colnames(xCell_result) %in% male_tumorFree]
    gene_female_tumor = as.numeric(xCell_result[index, colnames(xCell_result) %in% female_tumor])
    gene_female_tumorFree = as.numeric(xCell_result[index, colnames(xCell_result) %in% female_tumorFree])
    
    
    #input data for plot
    gene_data = data.frame(value=c(gene_male_tumorFree, gene_female_tumorFree, gene_male_tumor, gene_female_tumor), 
                           type=c(rep('male_tumorFree', length(gene_male_tumorFree)),rep('female_tumorFree', length(gene_female_tumorFree)), rep('male_tumor', length(gene_male_tumor)), rep('female_tumor', length(gene_female_tumor))))
    
    gene_data = mutate(gene_data, model=factor(type, levels=c('female_tumorFree', 'male_tumorFree','female_tumor', 'male_tumor')))
    
    #violin plot with box plot
    comparison=list(c("female_tumorFree", "male_tumorFree"),c("female_tumor", "male_tumor"))
    
    ggviolin(gene_data, x = "model", y = "value", fill = "model",
             #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "boxplot", add.params = list(fill = "white"))+
        stat_compare_means(aes(type), comparisons = comparison) + # Add significance levels
    #    +stat_compare_means(label.y = 3)        # Add global the p-value 
        ggtitle(gene)  + ylab("abundance")
    
    ggsave(paste0("./figure2/",as.character(gene),".png"), width=25, height = 15, units='cm', limitsize = FALSE)
    rm(comparison, gene, gene_data)
    rm(gene_male_tumor, gene_male_tumor_B, gene_male_tumor_C, gene_male_tumor_no, gene_male_tumorFree, gene_female_tumor, gene_female_tumorFree)
} #before running this code, check the folder file
rm(i, index, xCell_pval, xCell_result)
