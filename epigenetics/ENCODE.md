# ENCODE [(homepage)](https://www.encodeproject.org/)

# download ENCODE data
download query list in ENCODE homepage. [(example data)](https://www.encodeproject.org/experiments/ENCSR872WGW/)


example of query list file :  
``` 
"https://www.encodeproject.org/metadata/?type=Experiment&%40id=%2Fexperiments%2FENCSR872WGW%2F&option=raw"
https://www.encodeproject.org/files/ENCFF624DNH/@@download/ENCFF624DNH.fastq.gz
https://www.encodeproject.org/files/ENCFF157PAR/@@download/ENCFF157PAR.fastq.gz
https://www.encodeproject.org/files/ENCFF795EHY/@@download/ENCFF795EHY.fastq.gz
https://www.encodeproject.org/files/ENCFF121EPT/@@download/ENCFF121EPT.fastq.gz
``` 


```
xargs -L 1 curl -O -J -L < files.txt
```
