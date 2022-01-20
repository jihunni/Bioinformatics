# parsing a list (subsystem-gene) into a dataframe(gene-subsystem)
load('genesByHMRSubsystems.RData')

genesByHMRSubsystems_df = data.frame(matrix(nrow=0, ncol=2)) #initialization
colnames(genesByHMRSubsystems_df) = c("gene", "subsystem")
for(iteration in 1:length(genesByHMRSubsystems)[1]){
  print(paste0("iteration: ", iteration))
  split_list = unlist(genesByHMRSubsystems[iteration])
  #print(split_list)
  output_df = data.frame(gene=split_list, subsystem=rep(names(genesByHMRSubsystems[iteration]), length(split_list)))
  genesByHMRSubsystems_df = rbind(genesByHMRSubsystems_df, output_df)
}
rownames(genesByHMRSubsystems_df)=1:dim(genesByHMRSubsystems_df)[1]
save(genesByHMRSubsystems_df, file='genesByHMRSubsystems_df.RData')

# To map a gene into a pathway
pathway_score_df = data.frame(matrix(nrow=length(unique(gene_list_pathwayMatrix_2$subsystem)), ncol=2))
colnames(score_df) = c("subsystem", "score")
pathway_score_df$subsystem = unique(gene_list_pathwayMatrix_2$subsystem)
for(iteration in unique(gene_list_pathwayMatrix_2$subsystem)){
  subsystem_df = pathway_matrix_zscore[,colnames(pathway_matrix_zscore) %in% gene_list_pathwayMatrix_2$gene_symbol[gene_list_pathwayMatrix_2$subsystem==iteration]]
  score = sum(subsystem_df)/dim(subsystem_df)[2]
  if (length(score)>0){ #to escape the NA error
    pathway_score_df$score[score_df$subsystem==iteration] = score
  }
}
save(pathway_score_df, file='pathway_score_df_20220120.rda')
rm(score, subsystem_df, iteration)
