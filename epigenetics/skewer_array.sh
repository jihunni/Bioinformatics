#!/bin/bash
# run skewer as an array
# warning : add one vacant "" for input_output_array. The first element does not work properly.
# parameters
triming_seq="/opt/trimmomatic/0.39/adapters/NexteraPE-PE.fa"

function run_skewer {
	# function argument
	input_directory=$1
	output_directory=$2
	input_output_array=$3
	
	echo $input_directory
	for list_array in "${input_output_array[@]}" 
	do
		IFS=' ' read -r -a list <<< "$list_array"
		#list=$(echo ${list_array} | tr " " "\n")

		input1=${list[0]}
		input2=${list[1]}
		output_filename=${list[2]}
		
		#echo $input_output_array
		#echo ${input_output_array[@]}
		#echo $list_array
		#echo $input1
		#echo $input2
		#echo $output_filename		
		/opt/skewer/0.2.2/skewer-0.2.2-linux-x86_64 -f auto -t 15 -m pe -x ${triming_seq} -o ${output_directory}/${output_filename} $input_directory/${input1} $input_directory/${input2}
	done
}

#1
input_directory="${PWD}/01_ENCSR496PPU"
output_directory="$PWD/01_ENCSR496PPU/skewer"
input_output_array=(""  "ENCFF210LHI.fastq.gz ENCFF073GGX.fastq.gz ENCLB246GTQ_1" "ENCFF852QZF.fastq.gz ENCFF191IRU.fastq.gz ENCLB246GTQ_2" "ENCFF890VBZ.fastq.gz ENCFF284IIZ.fastq.gz ENCLB246GTQ_3" "ENCFF368WSR.fastq.gz ENCFF627ETL.fastq.gz  ENCLB246GTQ_4" "ENCFF064BBW.fastq.gz ENCFF899DAI.fastq.gz ENCLB875UFQ_1" "ENCFF579EPH.fastq.gz ENCFF075OLL.fastq.gz ENCLB875UFQ_2" "ENCFF827LXR.fastq.gz ENCFF391BZQ.fastq.gz ENCLB875UFQ_3" "ENCFF674VOA.fastq.gz ENCFF910VBI.fastq.gz ENCLB875UFQ_4")
run_skewer $input_directory $output_directory $input_output_array

#2
input_directory="${PWD}/02_ENCSR600ZHS"
output_directory="$PWD/02_ENCSR600ZHS/skewer"
input_output_array=("" "ENCFF448VTT.fastq.gz ENCFF010DJE.fastq.gz ENCLB347CEE_1" "ENCFF056PQW.fastq.gz ENCFF762DAM.fastq.gz ENCLB347CEE_2" "ENCFF435QMP.fastq.gz ENCFF082VIV.fastq.gz ENCLB347CEE_3" "ENCFF419FSN.fastq.gz ENCFF867KZO.fastq.gz ENCLB347CEE_4" "ENCFF250QRP.fastq.gz ENCFF965JGU.fastq.gz ENCLB555UOT_1" "ENCFF254GWH.fastq.gz ENCFF936PIY.fastq.gz ENCLB555UOT_2" "ENCFF470RAF.fastq.gz ENCFF540GNM.fastq.gz ENCLB555UOT_3" "ENCFF754MTD.fastq.gz ENCFF993LYP.fastq.gz ENCLB555UOT_4")
run_skewer $input_directory $output_directory $input_output_array
