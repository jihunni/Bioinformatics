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
	```

	```
