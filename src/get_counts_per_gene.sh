# Transform each BAM on classic BED (no use of -b)
bedtools bamtobed -i results/mapping/without_duplicates/2017_09_20_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.bam > results/mapping/without_duplicates/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_bam_to_bed.bed &&\
bash src/date.sh results/mapping/without_duplicates/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_bam_to_bed.bed &&\ 

bedtools bamtobed -i results/mapping/without_duplicates/2017_09_20_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.bam > results/mapping/without_duplicates/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_bam_to_bed.bed &&\
bash src/date.sh results/mapping/without_duplicates/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_bam_to_bed.bed &&\

# Make count of nb of reads falling within genes 
bedtools map -o count -a data/ReferenceGenomes/2017_12_06_Mesorhabditis_belari_JU2817_v2_genes_1based.bed -b results/mapping/without_duplicates/2017_11_12_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_bam_to_bed.bed > results/coverage/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_count_genes.txt &&\ 
bash src/date.sh results/coverage/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_count_genes.txt &&\

bedtools map -o count -a data/ReferenceGenomes/2017_12_06_Mesorhabditis_belari_JU2817_v2_genes_1based.bed -b results/mapping/without_duplicates/2017_11_12_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_bam_to_bed.bed > results/coverage/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_count_genes.txt &&\
bash src/date.sh results/coverage/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_count_genes.txt

