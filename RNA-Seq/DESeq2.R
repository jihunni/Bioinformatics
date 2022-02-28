library(tidyverse)
library(ggplot2)
library("ggpubr") #ggplot with p-value

load(propensity_normalized)
propensity_normalized = rbind(tumor_noHepatitis.propensity.data, tumorFree_noHepatitis.propensity.data)
propensity_normalized = propensity_normalized[propensity_normalized$weights != 0,]

HCC_protein_sample = rownames(HCC_protein)
HCC_protein2 = HCC_protein[,colnames(HCC_protein) %in% propensity_normalized$case]
rm(HCC_protein2)
HCC_protein =cbind(gene=HCC_protein_sample, HCC_protein)
colnames(HCC_protein)[1] = 'SampleID'

#statistics
statistics = propensity_normalized %>%
    group_by(sex, status) %>%
    summarise(number=n())

### Propensity Algorithm ###
weight_matrix = data.frame(case = colnames(HCC_protein)[-1])
weight_matrix = left_join(weight_matrix, propensity_normalized, by=c('case' = 'case'))
weight_matrix = select(weight_matrix, case, weights)
weight_vector = weight_matrix$weights
HCC_protein_propensity = data.frame(mapply('*',HCC_protein[,-1],weight_vector))
##HCC_protein_propensity[,2] == HCC_protein[,3] * weight_vector[2]
HCC_protein_propensity = cbind(data.frame(SampleID = HCC_protein[,1]), HCC_protein_propensity)
colnames(HCC_protein_propensity) = colnames(HCC_protein)
rownames(HCC_protein_propensity) = rownames(HCC_protein)
save(HCC_protein_propensity, file='HCC_protein_propensity_2021.05.08.rda')
rm(weight_vector, HCC_protein_sample, HCC_protein)

write.table(HCC_protein_propensity, file = "HCC_protein_propensity.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote=FALSE, append=FALSE)
#### protein level expression per gene ###
male_tumor = propensity_normalized$case[propensity_normalized$sex == 'Male' & propensity_normalized$status == 'With Tumor']
male_tumorFree = propensity_normalized$case[propensity_normalized$sex == 'Male' & propensity_normalized$status == 'Tumor Free']
female_tumor = propensity_normalized$case[propensity_normalized$sex == 'Female' & propensity_normalized$status == 'With Tumor']
female_tumorFree = propensity_normalized$case[propensity_normalized$sex == 'Female' & propensity_normalized$status == 'Tumor Free']

# Violin plots with box plots inside - for loop
for(i in HCC_protein_propensity$SampleID){
#    i = 'DJ1'
    print(i)
    gene = i

    #all group comparison between male and female
    gene_male_tumorFree = as.numeric(HCC_protein_propensity[gene, colnames(HCC_protein_propensity) %in% male_tumorFree])
    gene_male_tumor = as.numeric(HCC_protein_propensity[gene, colnames(HCC_protein_propensity) %in% male_tumor])
    gene_female_tumorFree = as.numeric(HCC_protein_propensity[gene, colnames(HCC_protein_propensity) %in% female_tumorFree])
    gene_female_tumor = as.numeric(HCC_protein_propensity[gene, colnames(HCC_protein_propensity) %in% female_tumor])
    
    #input data for plot
    gene_data = data.frame(value=c(gene_male_tumorFree, gene_male_tumor, gene_female_tumorFree, gene_female_tumor),
                           
                           type=c(rep('male_tumorFree', length(gene_male_tumorFree)), rep('male_tumor', length(gene_male_tumor)), rep('female_tumorFree', length(gene_female_tumorFree)), rep('female_tumor', length(gene_female_tumor))))
    
    gene_data = mutate(gene_data, model=factor(type, levels=c("female_tumorFree", 'male_tumorFree', 'female_tumor', 'male_tumor')))
    
    #violin plot with box plot
    comparison=list(c("female_tumorFree", "male_tumorFree"), c("female_tumor", "male_tumor"))
    ggviolin(gene_data, x = "model", y = "value", fill = "model",
             #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "boxplot", add.params = list(fill = "white"))+
        ggtitle(gene) + ylab("protein expression level") +
        stat_compare_means(aes(type), comparisons = comparison) # Add significance levels
    #    +stat_compare_means(label.y = 3)        # Add global the p-value 
    ggsave(paste0("./figure/",as.character(gene),".png"), width=25, height = 15, units='cm', limitsize = FALSE)
    rm(comparison, gene, gene_data, gene_female_tumor, gene_female_tumorFree, gene_male_tumor, gene_male_tumorFree)
} #before running this code, check the folder file
