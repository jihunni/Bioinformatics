# To prepare the download links for bacterial genome in RefSeq
# Filter : complete genome / reference genome, representative genome
library(readr)
RefSeq_list = read_tsv('RefSeq_bacterial_assembly_summary.txt', skip = 1)
colnames(RefSeq_list)[1]='assembly_accession'
table(RefSeq_list$refseq_category)
table(RefSeq_list$assembly_level)
#na      reference genome representative genome 
#222839                    15                 14994

RefSeq_list_complegeGenome = RefSeq_list[RefSeq_list$assembly_level == 'Complete Genome',]
table(RefSeq_list_complegeGenome$genome_rep)
table(RefSeq_list_complegeGenome$refseq_category)
#na      reference genome representative genome 
#21684                    ￣                  3582
table(is.na(RefSeq_list_complegeGenome$ftp_path))

output = RefSeq_list_complegeGenome[RefSeq_list_complegeGenome$refseq_category != 'na',]$ftp_path
output_fullAddress = lapply(output, function(x) {
  list = strsplit(x, '/')
  #print(unlist(list)[10])
  x = paste0(x,'/',unlist(list)[10],'_genomic.fna.gz')
})
output_fullAddress = unlist(output_fullAddress)
write.table(output_fullAddress, file = "RefSeq_bacterial_FTPpath.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote=FALSE, append=FALSE)
