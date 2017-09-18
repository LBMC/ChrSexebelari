# Pich header from BAM
samtools view -H results/mapping/without_duplicates/2017_09_13_MRDR5_trim_Mbelari_mapped_sort_rmdup.bam > tmp/2017_09_18_MRDR5_trim_Mbelari_mapped_sort_rmdup_header.sam &&\
samtools view -H results/mapping/without_duplicates/2017_09_13_MRDR6_trim_Mbelari_mapped_sort_rmdup.bam > tmp/2017_09_18_MRDR6_trim_Mbelari_mapped_sort_rmdup_header.sam &&\

# Remove lines within header which are not associated to M.belari contigs
awk -F "\t" '$1 == "@HD" || $2 ~ /^SN:M/ {print $0}' tmp/2017_09_18_MRDR5_trim_Mbelari_mapped_sort_rmdup_header.sam > tmp/2017_09_18_MRDR5_reheader.sam &&\
awk -F "\t" '$1 == "@HD" || $2 ~ /^SN:M/ {print $0}' tmp/2017_09_18_MRDR6_trim_Mbelari_mapped_sort_rmdup_header.sam > tmp/2017_09_18_MRDR6_reheader.sam &&\

# Change header
samtools reheader tmp/2017_09_18_MRDR5_reheader.sam results/mapping/without_duplicates/2017_09_13_MRDR5_trim_Mbelari_mapped_sort_rmdup.bam > results/mapping/without_duplicates/2017_09_18_MRDR5_trim_Mbelari_mapped_sort_rmdup.bam

samtools reheader tmp/2017_09_18_MRDR6_reheader.sam results/mapping/without_duplicates/2017_09_13_MRDR6_trim_Mbelari_mapped_sort_rmdup.bam > results/mapping/without_duplicates/2017_09_18_MRDR6_trim_Mbelari_mapped_sort_rmdup.bam



