module load bedops

sort-bed --tmpdir $PWD/tmp /home/jihun/data/ENCODE/HCT116_ENCSR872WGW/HMMRATAC/ENCLB310CQI_summits.bed > HCT116_ENCLB310CQI_summits.sorted.bed

sort-bed --tmpdir $PWD/tmp /home/jihun/data/ENCODE/HCT116_ENCSR872WGW/HMMRATAC/ENCLB404KGX_summits.bed > HCT116_ENCLB404KGX_summits.sorted.bed


module unload bedops
