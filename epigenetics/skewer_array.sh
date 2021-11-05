# run skewer as an array
## parameters
input_directory="${PWD}/12_ENCSR548QCP"
input_output_array=("ENCFF138HKP.fastq.gz ENCFF565VHZ.fastq.gz NCLB246GTQ") #(input1 input2 output_Name)
triming_seq="/opt/trimmomatic/0.39/adapters/NexteraPE-PE.fa"
output_directory="$PWD/01_ENCSR496PPU/skewer"

# run skewer
echo $input_directory
for list_array in "${input_output_array[@]}" 
do
	IFS=' ' read -r -a list <<< "$list_array"
	#list=$(echo ${list_array} | tr " " "\n")

	#echo "list0"
	#echo ${list[0]}
	#echo "list1"
	#echo ${list[1]}
	#echo "list2"
	#echo ${list[2]}
	input1=${list[0]}
	input2=${list[1]}
	output_filename=${list[2]}
		
	/opt/skewer/0.2.2/skewer-0.2.2-linux-x86_64 -f auto -t 10 -m pe -x ${triming_seq} -o ${output_directory}/${output_filename} $input_directory/${input1} $input_directory/${input2}
	
	#echo "1\n"
	#echo less ${output_directory}/${output_filename}
	#echo "2\n"
	#echo less $input_directory/${input1}
	#echo "3\n"
	#echo less $input_directory/${input2}
	 
done

