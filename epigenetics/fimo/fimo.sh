#!/bin/bash

#load module
module load MEME_Suite/5.4.1

echo "SLUM_ARRAY_TASK_ID : ${SLURM_ARRAY_TASK_ID}"

# execute FIMO
main () {
## set the directory
if (( ${SLURM_ARRAY_TASK_ID} < 10 )); then
	MOTIF_FILE=~/python/tfmotif_split/TRANSFAC_tfmatrix_human_MEME_00${SLURM_ARRAY_TASK_ID}.txt
	OUTPUT_DIRECTORY=~/data/motif/fimo_00${SLURM_ARRAY_TASK_ID}/
elif (( ${SLURM_ARRAY_TASK_ID} < 100 )); then
	MOTIF_FILE=~/python/tfmotif_split/TRANSFAC_tfmatrix_human_MEME_0${SLURM_ARRAY_TASK_ID}.txt
	OUTPUT_DIRECTORY=~/data/motif/fimo_0${SLURM_ARRAY_TASK_ID}/
else
	MOTIF_FILE=~/python/tfmotif_split/TRANSFAC_tfmatrix_human_MEME_${SLURM_ARRAY_TASK_ID}.txt
	OUTPUT_DIRECTORY=~/data/motif/fimo_${SLURM_ARRAY_TASK_ID}/
fi
SEQUENCE_FILE=/home/data/ref_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa

## run fimo
fimo  --oc ${OUTPUT_DIRECTORY} ${MOTIF_FILE} ${SEQUENCE_FILE}
}

time main

# unload module
module unload MEME_Suite/5.4.1
# declare the termination of a program
echo "finish"
