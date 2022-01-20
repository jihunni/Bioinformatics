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
