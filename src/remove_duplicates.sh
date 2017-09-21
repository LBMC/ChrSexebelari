# Remove duplicate
# Done locally due to problems with PSMN, using samtools 1.3

##### Args #####
# $1: input BAM
# $2: output BAM without PCR duplicates
# $3: flagstat file output on $2

##### remove of duplicates
samtools rmdup $1 $2 &&\

samtools flagstat $2 > $3
