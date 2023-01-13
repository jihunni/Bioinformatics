# merge all multiple output files from FIMO in a one tsv
# All output files are located in 'fimo_output' directory
output_fileName='merge_fimo_result.gff'
# iterate
cat header_file.gff >  ${output_fileName}

for file_name in $(ls fimo_output1|grep fimo_); do
	echo ${file_name} ;
	#cat fimo_output1/${file_name}/fimo.tsv >> merge_fimo_result_3.tsv 
		#problem : the first header line should be removed
	sed -n '2,$p' fimo_output1/${file_name}/fimo.gff >>  ${output_fileName}
done


for file_name in $(ls fimo_output2|grep fimo_); do
        echo ${file_name} ;
	sed -n '2,$p' fimo_output2/${file_name}/fimo.gff >>  ${output_fileName}
done

