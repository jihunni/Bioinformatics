# IO
input_file='merge_fimo_result.gff'
output_file_name='merge_fimo_result'

#load modules
module load bedops

#convert
## && : to run the second command only if the first exited successfully
gff2bed --max-mem 200G  < $input_file > ${output_file_name}.bed && sort-bed --max-mem 200G  ${output_file_name}.bed > ${output_file_name}.sorted.bed
