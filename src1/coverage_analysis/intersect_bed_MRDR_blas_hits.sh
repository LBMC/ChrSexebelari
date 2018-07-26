awk 'OFS="\t" { print $1, $2, $2, $3 }' ~/Documents/stage_mbelari/results/coverage_analysis/2018_06_21_MRDR5_vs_Hybrid_assembly.sorted.filtered.bed > tmp_input_intersect_d_genes.bed &&\
bedtools intersect -a tmp_input_intersect_d_genes.bed -b ~/Documents/stage_mbelari/results/annotation/coverage_analysis/Mesorhabditis_belari_JU2817_hybrid_assembly_genes_1based.bed -wb > ~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-17-MRDR5_vs_Hybrid_assembly_sort_in_genes_bp.bed 

awk 'OFS="\t" { print $1, $2, $2, $3 }' ~/Documents/stage_mbelari/results/coverage_analysis/2018_06_21_MRDR6_vs_Hybrid_assembly.sorted.filtered.bed > tmp_input_intersect_d_genes.bed &&\
bedtools intersect -a tmp_input_intersect_d_genes.bed -b ~/Documents/stage_mbelari/results/annotation/coverage_analysis/Mesorhabditis_belari_JU2817_hybrid_assembly_genes_1based.bed -wb > ~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-17-MRDR6_vs_Hybrid_assembly_sort_in_genes_bp.bed 
