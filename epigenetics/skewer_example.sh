# shell script runned by jihun_workstation
# example

input_directory="${PWD}/01_ENCSR496PPU/"
input_1="ENCLB246GTQ_total_forward_paired.fastq"
input_2="ENCLB246GTQ_total_reverse_paired.fastq"

output_directory="$PWD/01_ENCSR496PPU/skewer/"
output_filename="ENCLB246GTQ"

/opt/skewer/0.2.2-linux-x86_64/skewer-0.2.2-linux-x86_64 -f auto -t 12 -m pe -x AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -o ${output_directory}/${output_filename} $input_directory/${input_1} $input_directory/${input_2} 
