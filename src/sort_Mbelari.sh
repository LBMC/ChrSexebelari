# Sort and index BAM files
# Done locally due to problems with PSMN, using samtools 1.3

##### Args #####
# $1: input BAM
# $2: output BAM 
# $3: output report with idxstat

samtools sort $1 -o $2 &&\ 
samtools index $2 &&\

##### idxstats
samtools idxstats $2 > $3
