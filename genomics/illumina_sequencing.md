# Adapter trimming: Why are adapter sequences trimmed from only the 3' ends of reads?
Ref : https://support.illumina.com/bulletins/2016/04/adapter-trimming-why-are-adapter-sequences-trimmed-from-only-the--ends-of-reads.html
09/03/21

Illumina FASTQ file generation pipelines include an adapter trimming option for the removal of adapter sequences from the 3’ ends of reads. Adapter sequences should be removed from reads because they interfere with downstream analyses, such as alignment of reads to a reference. The adapters contain the sequencing primer binding sites, the index sequences, and the sites that allow library fragments to attach to the flow cell lawn. Libraries prepared with Illumina library prep kits require adapter trimming only on the 3’ ends of reads, because adapter sequences are not found on the 5’ ends.

Note: Libraries prepared with the Nextera™ Mate Pair library prep kit are an exception, and guidelines for trimming adapters from these libraries can be found in the Data Processing of Nextera Mate Pair Reads on Illumina Sequencing Platforms technical note.

To understand why adapter sequences are found only on the 3’ ends of the reads, it helps first to understand where the sequencing primers anneal to the library template on a flow cell. The diagrams below show the sites of primer annealing at each stage of sequencing run: Read 1, Index 1, Index 2 and Read 2.

![image](https://user-images.githubusercontent.com/48517782/140247532-ffa951ad-66f5-4096-97b6-793f996ed4a6.png)  
Figure 1. MiSeq™, HiSeq™ 1000/1500/2000/2500 and NovaSeq™ 6000 v1.0 reagents paired-end flow cell

![image](https://user-images.githubusercontent.com/48517782/140247549-c25e48f8-a115-4a07-a27c-0e8caa89083c.png)  
Figure 2. iSeq™ 100, MiniSeq™, NextSeq™ 500/550, NextSeq 1000/2000, HiSeq 3000/4000 paired-end flow cell, and NovaSeq 6000 v1.5 reagents

As shown in Figures 1 and 2, in both Read 1 and Read 2, the sequencing primer anneals to the adapter, immediately upstream of the DNA insert (in gray). Because the sequencing starts at the first base of the DNA insert in Reads 1 and 2, the adapter is not sequenced at the start of the read. However, if the sequencing extends beyond the length of the DNA insert, and into the adapter on the opposite end of the library fragment, that adapter sequence will be found on the 3’ end of the read. Therefore, reads require adapter trimming only on their 3’ ends.
