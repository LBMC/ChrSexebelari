
# Transform each BAM on classic BED (no use of -b)
#bedtools bamtobed -i ~/Documents/stage_mbelari/results/13genes/2018-06-12-MRDR5_vs_Hybrid_assembly.sorted.bam  > 2018-07-13-MRDR5_vs_Hybrid_assembly.sorted.bam_to_bed.bed

#bedtools bamtobed -i ~/Documents/stage_mbelari/results/13genes/2018-06-12-MRDR6_vs_Hybrid_assembly.sorted.bam  > 2018-07-13-MRDR6_vs_Hybrid_assembly.sorted.bam_to_bed.bed

# Make count of nb of reads falling within genes 
#bedtools map -o count -a ~/Documents/stage_mbelari/results/annotation/BRAKER/soft_masked/mbelari_hybrid_cds_1.bed -b 2018-07-13-MRDR5_vs_Hybrid_assembly.sorted.filtered.bam_to_bed.bed > MRDR5_vs_Hybrid_assembly_count_genes.txt 

#bedtools map -o count -a ~/Documents/stage_mbelari/results/annotation/BRAKER/soft_masked/mbelari_hybrid_cds_1.bed -b 2018-07-13-MRDR6_vs_Hybrid_assembly.sorted.filtered.bam_to_bed.bed > MRDR6_vs_Hybrid_assembly_count_genes.txt 


#bedtools sort -i 2018-07-13-MRDR5_vs_Hybrid_assembly.sorted.bam_to_bed.bed -faidx ~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt  > 2018-07-13-MRDR5_vs_Hybrid_assembly.sorted.bam_to_bed_idx.bed

#bedtools sort -i 2018-07-13-MRDR6_vs_Hybrid_assembly.sorted.bam_to_bed.bed -faidx ~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt  > 2018-07-13-MRDR6_vs_Hybrid_assembly.sorted.bam_to_bed_idx.bed

bedtools intersect -a ~/Documents/stage_mbelari/results/annotation/BRAKER/soft_masked/mbelari_hybrid_cds_1_sorted.bed -b 2018-07-13-MRDR5_vs_Hybrid_assembly.sorted.bam_to_bed.bed -c -sorted > MRDR5_vs_Hybrid_assembly_count_genes.txt 

bedtools intersect -a ~/Documents/stage_mbelari/results/annotation/BRAKER/soft_masked/mbelari_hybrid_cds_1_sorted.bed -b 2018-07-13-MRDR6_vs_Hybrid_assembly.sorted.bam_to_bed.bed -c -sorted > MRDR6_vs_Hybrid_assembly_count_genes.txt 
