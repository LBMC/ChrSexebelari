# Get bed only for gene regions 
awk '{ if ($8 == "gene") {print $0}}' data/ReferenceGenomes/2017_10_26_Mesorhabditis_belari_JU2817_v2.bed > data/ReferenceGenomes/2017_11_12_Mesorhabditis_belari_JU2817_v2_genes.bed &&\

# Transform each BAM on classic BED (no use of -b)
bedtools bamtobed -i results/mapping/without_duplicates/2017_09_20_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.bam > results/mapping/without_duplicates/2017_11_12_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_bam_to_bed.bed &&\

bedtools bamtobed -i results/mapping/without_duplicates/2017_09_20_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.bam > results/mapping/without_duplicates/2017_11_12_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_bam_to_bed.bed &&\

# Make count of nb of reads falling within genes 
bedtools map -o count -a data/ReferenceGenomes/2017_11_12_Mesorhabditis_belari_JU2817_v2_genes.bed -b results/mapping/without_duplicates/2017_11_12_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_bam_to_bed.bed > results/coverage/2017_11_12_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_count_genes.txt &&\ 

bedtools map -o count -a data/ReferenceGenomes/2017_11_12_Mesorhabditis_belari_JU2817_v2_genes.bed -b results/mapping/without_duplicates/2017_11_12_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_bam_to_bed.bed > results/coverage/2017_11_12_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_count_genes.txt 

