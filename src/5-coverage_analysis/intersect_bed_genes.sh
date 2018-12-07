# From Claire Analysis
# This scripts makes an intersection of the entire bed coverage file at bp level on the genic region
# Run on local computer

#awk 'OFS="\t" { print $1, $2, $2, $3 }' results/annotation/gene_counts/2018-07-13-MRDR5_vs_Hybrid_assembly.sorted.bam_to_bed.bed > tmp_input_intersect_d_genes.bed &&\
#bedtools intersect -a tmp_input_intersect_d_genes.bed -b results/annotation/coverage_analysis/Mesorhabditis_belari_JU2817_hybrid_assembly_genes_1based.bed -wb > results/annotation/coverage_analysis/2018-07-16-MRDR5_vs_Hybrid_assembly_sort_in_genes.bed &&\

#awk 'OFS="\t" { print $1, $2, $2, $3 }' results/annotation/gene_counts/2018-07-13-MRDR6_vs_Hybrid_assembly.sorted.bam_to_bed.bed > tmp_input_intersect_d_genes.bed &&\
#bedtools intersect -a tmp_input_intersect_d_genes.bed -b results/annotation/coverage_analysis/Mesorhabditis_belari_JU2817_hybrid_assembly_genes_1based.bed -wb > results/annotation/coverage_analysis/2018-07-16-MRDR6_vs_Hybrid_assembly_sort_in_genes.bed

awk 'OFS="\t" { print $1, $2, $2, $3 }' results/coverage_analysis/2018_07_25_MRDR5_vs_Hybrid_assembly.sorted.filtered.bed > tmp_input_intersect_d_genes.bed
bedtools intersect -a tmp_input_intersect_d_genes.bed -b results/coverage_analysis/Mesorhabditis_belari_JU2817_hybrid_assembly_genes_1based.bed -wb > results/coverage_analysis/2018-07-25-MRDR5_vs_Hybrid_assembly_sort_in_genes_bp.bed

awk 'OFS="\t" { print $1, $2, $2, $3 }' results/coverage_analysis/2018_07_25_MRDR6_vs_Hybrid_assembly.sorted.filtered.bed > tmp_input_intersect_d_genes.bed
bedtools intersect -a tmp_input_intersect_d_genes.bed -b results/coverage_analysis/Mesorhabditis_belari_JU2817_hybrid_assembly_genes_1based.bed -wb > results/coverage_analysis/2018-07-25-MRDR6_vs_Hybrid_assembly_sort_in_genes_bp.bed
