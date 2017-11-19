# bedtools 2.22 on local
awk 'OFS="\t" { print $1, $2, $2, $3 }' results/coverage/2017_10_26_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_sort.bed > results/coverage/tmp_input_intersect_d_genes.bed &&\
bedtools intersect -a results/coverage/tmp_input_intersect_d_genes.bed -b data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_genes.bed -wb > results/coverage/2017_11_18_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_sort_in_genes.bed &&\

awk 'OFS="\t" { print $1, $2, $2, $3 }' results/coverage/2017_10_26_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_sort.bed > results/coverage/tmp_input_intersect_d_genes.bed &&\
bedtools intersect -a results/coverage/tmp_input_intersect_d_genes.bed -b data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_genes.bed -wb > results/coverage/2017_11_18_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_sort_in_genes.bed 

