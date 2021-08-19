# sequence alignment overview
- procedure   
  <img src="https://user-images.githubusercontent.com/48517782/129913198-5ac027d0-082b-4ab9-b2a8-2c570af08b3f.png" width="200px" height="300px" title="" alt=""></img><br/>
- objective  
  ![image](https://user-images.githubusercontent.com/48517782/129914645-c7332904-65a3-425b-94fa-7243735532ec.png){: width="300" height="300"}{: .center}


# Data preparation
- NCBI Reference genome | [Mus musculus (house mouse)](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27)   
  ![image](https://user-images.githubusercontent.com/48517782/129909944-7ba6c7fb-c94a-4e1a-8351-9d81ad775069.png)   

- ensembl reference genome | [Mus musculus (house mouse)](https://asia.ensembl.org/Mus_musculus/Info/Index) | [ftp](http://ftp.ensembl.org/pub/release-103/fasta/mus_musculus/dna/)   
  ![image](https://user-images.githubusercontent.com/48517782/129910326-c302859e-f84b-4ebc-88e2-4f89f94785d2.png)
  file:
    - Mus_musculus.GRCm39.dna.toplevel.fa.gz 
    - Mus_musculus.GRCm39.dna_rm.toplevel.fa.gz 
    - Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz  

      > Unmasked : general purpose   
      > Masked : loss of information  
      > Soft-masked  
      > [reference](https://genestack.com/blog/2016/07/12/choosing-a-reference-genome/)  


  ```
  nohup wget http://ftp.ensembl.org/pub/release-103/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz 2> wget_reference_genome.out &
  ```

- sample
```
$fasterq-dump --split-files SRR11180057

$
for ((i=30 ; i<=80 ; i++ ))
do
	../../program/sratoolkit.2.11.0-ubuntu64/bin/fastq-dump --split-files SRR9141$i	
	echo 'done $i'
done
```

# fastQC : quality check
- install fastqc [reference](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt)
```
$ chmod 755 fastqc
$ ./fastqc
$ sudo ln -s /path/to/FastQC/fastqc /usr/local/bin/fastqc
```
- linux server supports graphic version
![image](https://user-images.githubusercontent.com/48517782/129911884-3501a837-0bde-44fa-a43d-85ca247bbd15.png)
![image](https://user-images.githubusercontent.com/48517782/129912021-e18119ee-7a31-40e4-90f6-92537a25e55c.png)

# cudadapt [link](https://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming)
```
$ cutadapt -q 10 -o output.fastq input.fastq
```

# STAR alignment


```
$ /home/bioinfo20165164/program/STAR-2.7.8a/source/STAR --runThreadN 8 \
	--quantMode TranscriptomeSAM GeneCounts \
	--genomeDir /home/bioinfo20165164/data/reference/Mus_musculus.GRCm39.dna.toplevel/star_index/ \
	--readFilesIn /home/bioinfo20165164/data/T-cell/raw_data/SRR914062.fastq \
	--outFileNamePrefix /home/bioinfo20165164/data/T-cell/result_STAR_alignment\star_SRR914062_ \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattributes Standard
```
> options:
>* --runThreadN: number of threads  
>* --runMode: genomeGenerate mode  
>* --genomeDir: /path/to/store/genome_indices  
>* --genomeFastaFiles: /path/to/FASTA_file  
>* --sjdbGTFfile: /path/to/GTF_file  
>* --sjdbOverhang: readlength -1  

## samtools index
- Indexing a genome sorted BAM file allows one to quickly extract alignments overlapping particular genomic regions. Moreover, indexing is required by genome viewers such as IGV so that the viewers can quickly display alignments in each genomic region to which you navigate.
- If an output filename is given, the index file will be written to out.index.

```
$ samtools index sample.sorted.bam
```
```
$ /home/bioinfo20165164/program/samtools-1.12/samtools index -@ 6 /home/bioinfo20165164/data/T-cell/result_STAR_alignment/GSM1169465_Aligned.sortedByCoord.out.bam

$for ((i=91 ; i<=99 ; i++))
do
	/home/bioinfo20165164/program/samtools-1.12/samtools index -@ 6 ./result_STAR_alignment/GSM11694${i}_Aligned.sortedByCoord.out.bam
	echo done ${i}
done
```

## samtools sort
```
$ samtools sort star_.bam -o SRR.sam -@ 6
```
  > options:
  > -o output file
  > -@ int #set number of sorting and compression threads

  output example: 
  ![image](https://user-images.githubusercontent.com/48517782/129912698-4dfa6d46-7ea2-43d7-b798-ca82983fde7e.png)

# HTSeq : to convert BAM into a gene count tsv file
install
https://htseq.readthedocs.io/en/master/install.html

HTSeq count command: 
```
$ python -m HTSeq.scripts.count [options] <alignment_files> <gff_file>
```
Output example:
![image](https://user-images.githubusercontent.com/48517782/129913899-9194a047-bd5b-48aa-a092-319ba818c767.png)
![image](https://user-images.githubusercontent.com/48517782/129913932-57d2db3c-36a8-4225-b6b4-81e10bbcb3b1.png)



# kallisto
```
$ kallisto index -i homo_sapiens_index /group/bioinforamtics.2021/reference/homo_sapiens/Homo_sapiens.GRCh38.cdna.all.fa
```
```
$  kallisto quant -i homo_sapiens_index -o output_kallisto quiz_1.cutadapt.fastq quiz_2.cutadapt.fa 
```

# Reference
- https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html