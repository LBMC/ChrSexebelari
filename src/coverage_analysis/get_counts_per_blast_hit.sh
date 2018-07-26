#From claire analysis

# Obtain raw counts per gene region: this require to transform BAM to BED and to intersect the BED with the annotation of the genic regions in BED format
# Run on local computer

# Transform each BAM on classic BED (no use of -b)
#bedtools bamtobed -i results/mapping/without_duplicates/2017_09_20_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.bam > results/mapping/without_duplicates/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_bam_to_bed.bed &&\
#bash src/date.sh results/mapping/without_duplicates/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_bam_to_bed.bed &&\ 

#bedtools bamtobed -i results/mapping/without_duplicates/2017_09_20_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.bam > results/mapping/without_duplicates/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_bam_to_bed.bed &&\
#bash src/date.sh results/mapping/without_duplicates/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_bam_to_bed.bed &&\

#Sort bed
#-faidx ~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt

#bedtools sort -faidx ~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt -i ~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-24-Blast_female_absent_genes_hits.bed > ~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-24-Blast_female_absent_genes_hits.sorted.bed

# Make count of nb of reads falling within genes -g ~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt
bedtools map -o count -g ~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt -a ~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-24-Blast_female_absent_genes_hits.sorted.bed -b ~/Documents/stage_mbelari/results/coverage_analysis/2018-06-12-MRDR5_vs_Hybrid_assembly.sorted.filtered.bam_to_bed.bed > ~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-24-Blast_female_absent_genes_hits_MRDR5_counts.txt &

bedtools map -o count -g ~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt -a ~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-24-Blast_female_absent_genes_hits.sorted.bed -b ~/Documents/stage_mbelari/results/coverage_analysis/2018-06-12-MRDR6_vs_Hybrid_assembly.sorted.filtered.bam_to_bed.bed > ~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-24-Blast_female_absent_genes_hits_MRDR6_counts.txt &

