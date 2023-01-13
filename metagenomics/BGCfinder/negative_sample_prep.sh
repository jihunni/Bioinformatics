# a program that prepares the negative set for bacterial genome mining.
REF_GENOME_DIR='/home/data/ref_genome/bacterial_genome/'
REF_GENOME_EXTENSION='fna'
BLAST_DB_DIR='/home/jihun/data/blastDB_mibig2.0_prot/mibig_prot_seqs_2.0.fasta'
OUTPUT_DIR='/home/jihun/data/bacterial_protein_negative/output/'
NUM_THREADS=25

module load prodigal/2.6.3
module load blast

#module load anaconda
#conda init bash
#conda activate bioPython

for file_name in $(ls ${REF_GENOME_DIR}); do
	file_name="${file_name##*/}" # to get fileName with extension (to remove path)
	file_name="${file_name%.*}" # to remove file extenstion
	echo ${file_name};
	
	# to continue the previous iteration
	if [ -e ${OUTPUT_DIR}/${file_name}.blast_6 ] ;
		then
		echo 'blast output result file exists'
	else
		
		# gene prediction
		prodigal -i ${REF_GENOME_DIR}/${file_name}.${REF_GENOME_EXTENSION} -a ${OUTPUT_DIR}/${file_name}.faa -o prodigal.output
		echo 'log::prodigal is done'	

		echo 'log::blastp starts'
		# To prepare the filtering : 1. blast
		blastp -db ${BLAST_DB_DIR} -query ${OUTPUT_DIR}/${file_name}.faa -outfmt 6 -num_threads ${NUM_THREADS} -out ${OUTPUT_DIR}/${file_name}.blast_6
		echo 'log::blastp is done'
	fi
	
	# Filtering
	python filter_reads_in_blast.py --input-fasta-file ${OUTPUT_DIR}/${file_name}.faa --blast-result-file ${OUTPUT_DIR}/${file_name}.blast_6 --output-fasta-file ${OUTPUT_DIR}/${file_name}.faa_filtered
	
	#break # to test
done
echo 'log::a program is finished'
