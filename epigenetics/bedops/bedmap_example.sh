# load module
module load bedops

#input
reference_file='reference_bedmap_motifs_2.bed'
data_file='reference_bedmap_map_2.bed'
#data_file='/home/jihun/data/TCGA-COAD_ATAC/HMMRATAC_slurm_2/eedca752-bfb9-4a73-a800-613aa1f5fe27_atacseq_gdc_realn_summits.sorted.bed'

#output
output_file='bedmap_example.bed'
# --echo-map flag gathers overlapping mapped elements for every reference elements
bedmap --echo  --mean ${reference_file} ${data_file} > $output_file
#bedmap --chrom chr21 --skip-unmapped --echo-map ${reference_file} ${data_file} > example_out.bed
