#From claire analysis

# Obtain raw counts per gene region: this require to transform BAM to BED and to intersect the BED with the annotation of the genic regions in BED format
# Run on local computer

# Transform each BAM on classic BED (no use of -b)
bedtools bamtobed -i results/coverage_analysis/2018-06-12-MRDR5_vs_Hybrid_assembly.sorted.filtered.bam > \
  results/2018-06-12-MRDR5_vs_Hybrid_assembly.sorted.filtered.bam_to_bed.bed
#bash src/date.sh results/mapping/without_duplicates/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_bam_to_bed.bed &&\

bedtools bamtobed -i results/coverage_analysis/2018-06-12-MRDR6_vs_Hybrid_assembly.sorted.filtered.bam > \
  results/2018-06-12-MRDR6_vs_Hybrid_assembly.sorted.filtered.bam_to_bed.bed
#bash src/date.sh results/mapping/without_duplicates/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_bam_to_bed.bed &&\

#Sort bed
bedtools sort -faidx results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt -i results/annotation/Mesorhabditis_belari_JU2817_hybrid_assembly_genes_1based.bed > results/annotation/Mesorhabditis_belari_JU2817_hybrid_assembly_genes_1based.sorted.bed

# Make count of nb of reads falling within genes
bedtools map -o count \
  -g results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt \
  -a results/annotation/Mesorhabditis_belari_JU2817_hybrid_assembly_genes_1based.sorted.bed \
  -b results/2018-06-12-MRDR5_vs_Hybrid_assembly.sorted.filtered.bam_to_bed.bed > \
  results/coverage_analysis/2018-07-24-MRDR5_vs_Hybrid_assembly.sorted.count.genes.txt

bedtools map -o count \
  -g results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt \
  -a results/annotation/Mesorhabditis_belari_JU2817_hybrid_assembly_genes_1based.sorted.bed \
  -b results/2018-06-12-MRDR6_vs_Hybrid_assembly.sorted.filtered.bam_to_bed.bed > \
  results/coverage_analysis/2018-07-24-MRDR6_vs_Hybrid_assembly.sorted.count.genes.txt



