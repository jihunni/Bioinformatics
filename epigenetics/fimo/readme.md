# FIMO
## Input
- transfac2meme_human_6715_oneExample.txt
- reference genome

## To run `FIMO` with multiple cores
```
fimo_slurm_array.sh
  fimo.sh
```

## output
- cisml.xml
- fimo.gff
- fimo.html
- fimo.tsv : fimo_output_tsv_example.tsv
- fimo.xml

## To merge all fimo output files
```
convert_gff_to_bed.sh   # to convert gff file to bed file
merge_fimo_tsv_result.sh
  header.tsv
merge_fimo_xml_result.sh
merge_fimo_gff_result.sh
  header_file.gff
```
gff -[manual conversion]-> bed  -[bedops]-> processed bed files
