# metagenomics
# Database
- JGI IMG : metagenomics
- GreenGenes : 16S
- SILVA
- RDP
## reference genome
- Genbank : 16S
- NCBI RefSeq : reference genome
	- To download reference bacterial genome, 
		- https://github.com/jihunni/Bioinformatics/blob/main/metagenomics/RefSeq_ftp_download_prep.R
		- https://github.com/jihunni/Bioinformatics/blob/main/metagenomics/RefSeq_ftp_download.sh
- EMBL-EBI
	ref : https://www.ebi.ac.uk/genomes/bacteria.html
	
## secondary metabolite
- MIBiG : BiG
- antiSMASH : BiG

## Genebank
1) Using the batch Entrez website  
https://www.ncbi.nlm.nih.gov/sites/batchentrez  
2) Using Perl: (copy into your terminal and press return/enter)
	```
	perl -e 'use LWP::Simple;getstore("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=".join(",",qw(6701965 6701969 6702094 6702105 6702160)),"seqs.fasta");'
	```
# Microbial community analysis
## Quality Control of raw reads
### Quality Control
- multiQC
### trimming
- Trimmomatic
  - fast
  - paired end data
- PRINSEQ
- FastX
  - paired-end data X
  - simple
- TagCleaner

## Combine paired reads to contigs
-  Mothur
	Ref : https://mothur.org/   
	-  Input : Tar package of FASTQ files
	-  Creates a reverse complement of the reverse read
	-  Performs a Neeleman alignment for the two reads
	- output : 
		- contigs.fasta.gz
		- samples.fastqs.txt
		- contigs.groups
		- contig.numbers.txt
		- contigs.summary.tsv

## Filter contigs and remove identical sequences
- Mothur `screens.seqs` : to filter contigs
	- length, ambigious bases, homopolymer
	- Input : contigs fasta file, groups file
	- Output : 
		- screened.fasta.gz
		- screened.groups
		- summary.screened.tsv
- Extract Unique sequences : to remove identical sequences
	- Input : fasta file, group files
	- Output : 
		- unique.fasta
		- unique.count_table
		- unique.summary.tsv
## Align sequences to reference template alginment
- SILVA reference template alginment
	Ref : https://mothur.org/wiki/silva_reference_files/
	- Input : unique.fasta.gz, unique.count_table
	- Steps:
		- find the cloest template sequence for the query sequence using K-meer search with 8mers
		- align the query and the de-gapped template sequence using Needleman-Wunsch pairwise alignment
		- re-insert gaps to the query and template pairwise alignment using the NAST algorithmn so that the query sequence alignment is compatible with the original template alignment
	- Ouptut:
		- aligned.fasta.gz
		- custom.reference.summary.tsv
		- aligned-summary.tsv  

## Filter and trim aligned sequence
- filter
- trim sequence alignment for overhangs and empty columns
	- `Filter sequence alignment`
	- Input : screened.fasta.gz, screend.count_table
	- Output :
		- filtered-unqiue.fasta.gz
		- filtered-unique.count_table
		- filtered-unique-summary.tsv
		- filtered-log.txt  

## Remove sequencing errors and chimeras
- `Precluster algined sequences`
	- Input : filtered-unique.fasta.gz, filtered-unique.count_table
	- Output : 
- Chimeras : artifact sequence formed by two biological sequences. (e.g. incomplete extension during PCR)
	- `Remove chimeric sequences`
	-	Input : 
	-	Output : chimeras.removed.fasta.gz, chimeras.removed.count_table, chimeras.removed.summary.tsv

## Classify sequences to taxonomic units
- Silva, UNITE
- Wang method : probabilistic method
- Input : chimeras.removed.fasta.gz, chimeras.removed.count_table
- Output : sequences-taxonomy-assignment.txt, classification-summary.tsv
- to remove unwanted lineages (e.g. MT, chloroplasts, unknown)
