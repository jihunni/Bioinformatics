# conda activate bioPython

import argparse
import os
from Bio import SeqIO
from Bio import SearchIO

parser = argparse.ArgumentParser(description='')
parser.add_argument('--input-fasta-file', '-i', required=True, type=str, default=None, help='Input amino acid fasta file')
parser.add_argument('--blast-result-file', '-b', required=True, type=str, default=None, help='blast result file with format 6')
parser.add_argument('--output-fasta-file', '-o',  required=True, type=str, default=None, help='output fasta (amino acid) file')

args = parser.parse_args() #initialize
INPUT_FASTA_FILE = args.input_fasta_file # input amino acid fasta file
BLAST_RESULT = args.blast_result_file # result from blast (form 6)
OUTPUT_FASTA_FILE = args.output_fasta_file # the name assigned for output fasta amino acid file

# To make a ID list that are included in blast reads 
# (These reads are matched with any read in blastDB)
blast_records = SearchIO.parse(BLAST_RESULT, 'blast-tab')
blast_id_list = []
for blast_record in blast_records:
    blast_id_list.append(blast_record.id)
    
# To filter from a ID list that are included in blast reads 
# (These filtered reads are not matched with any read in blastDB)
seq_records = SeqIO.parse(INPUT_FASTA_FILE,"fasta")
filtered_sequences = []  # Setup an empty list
for record in seq_records:
    if not (record.id in blast_id_list):
        # If record is not matched with blast, add this record to our list
        filtered_sequences.append(record)
print("Found %i short sequences" % len(filtered_sequences))

# to write a file
with open(OUTPUT_FASTA_FILE, "w") as file:
    SeqIO.write(filtered_sequences, file, "fasta")
print("done")
