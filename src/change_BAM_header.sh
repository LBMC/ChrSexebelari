# Pich header from BAM after having added RG information and save it under SAM extension
mkdir tmp &&\
samtools view -H results/mapping/without_duplicates/2017_09_20_MRDR5_trim_Mbelari_mapped_rmdup_rg.bam > tmp/MRDR5_trim_Mbelari_mapped_rmdup_rg_header.sam &&\
samtools view -H results/mapping/without_duplicates/2017_09_20_MRDR6_trim_Mbelari_mapped_rmdup_rg.bam > tmp/MRDR6_trim_Mbelari_mapped_rmdup_rg_header.sam &&\

# Remove lines within header which are not associated to M.belari contigs
awk -F "\t" '$1 == "@HD" || $1 == "@RG" || $2 ~ /^SN:M/ {print $0}' tmp/MRDR6_trim_Mbelari_mapped_rmdup_rg_header.sam > results/mapping/without_duplicates/MRDR6_reheader.sam &&\
awk -F "\t" '$1 == "@HD" || $1 == "@RG" || $2 ~ /^SN:M/ {print $0}' tmp/MRDR5_trim_Mbelari_mapped_rmdup_rg_header.sam > results/mapping/without_duplicates/MRDR5_reheader.sam &&\

# Change header
samtools reheader results/mapping/without_duplicates/MRDR5_reheader.sam results/mapping/without_duplicates/2017_09_20_MRDR5_trim_Mbelari_mapped_rmdup_rg.bam > results/mapping/without_duplicates/2017_09_20_MRDR5_trim_Mbelari_mapped_rmdup_rgh.bam &&\
samtools reheader results/mapping/without_duplicates/MRDR6_reheader.sam results/mapping/without_duplicates/2017_09_20_MRDR6_trim_Mbelari_mapped_rmdup_rg.bam > results/mapping/without_duplicates/2017_09_20_MRDR6_trim_Mbelari_mapped_rmdup_rgh.bam

