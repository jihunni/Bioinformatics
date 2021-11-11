# Build a index file
Ref: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer
Ref: https://m.blog.naver.com/PostView.naver?isHttpsRedirect=true&blogId=ruleoutlife&logNo=220743247716

- general code
  ```
  bowtie2-build --threads [thread_number] [ref_genome.fasta] [output_name]
  ```

- running code
  ```
  # load module
  module load bowtie2/2.4.4

  # indexing
  bowtie2-build --threads 10 Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38
  ```

# alignment
	Notice that index file directory sholud be assigned on an environment variable `BOWTIE2_INDEXES`  
	```
	export BOWTIE2_INDEXES='/home/data/ref_genome/'
	bowtie2 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -X 2000 --rg-id -p 10 -x Homo_sapiens.GRCh38_bowtie2  -1 NCLB246GTQ-trimmed-pair1.fastq -2 NCLB246GTQ-trimmed-pair2.fastq -S output.sam
	```
	> option :  
 	> - p : the number of threads
	> - x : index file name (prefix)
	> -1 and -2 : paired-end input fastq file
