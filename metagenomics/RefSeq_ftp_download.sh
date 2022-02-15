#!/bin/bash
# download fasta file from a input list
fileName=/home/data/ref_genome/bacterial_genome/RefSeq_bacterial_FTPpath.txt
num=1
while read line; do
        echo "$line"
        echo $num
        wget ${line}
        let num++
done < $fileName
