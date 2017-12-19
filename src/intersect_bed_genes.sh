# bedtools 2.22 on local
# This scripts makes an intersection of the entire bed coverage file at bp level on the genic region

awk 'OFS="\t" { print $1, $2, $2, $3 }' results/coverage/2017_10_26_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_sort.bed > results/coverage/tmp_input_intersect_d_genes.bed &&\
bedtools intersect -a results/coverage/tmp_input_intersect_d_genes.bed -b data/ReferenceGenomes/2017_12_06_Mesorhabditis_belari_JU2817_v2_genes_1based.bed -wb > results/coverage/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_sort_in_genes.bed &&\
bash src/date.sh results/coverage/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_sort_in_genes.bed &&\

awk 'OFS="\t" { print $1, $2, $2, $3 }' results/coverage/2017_10_26_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_sort.bed > results/coverage/tmp_input_intersect_d_genes.bed &&\
bedtools intersect -a results/coverage/tmp_input_intersect_d_genes.bed -b data/ReferenceGenomes/2017_12_06_Mesorhabditis_belari_JU2817_v2_genes_1based.bed -wb > results/coverage/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_sort_in_genes.bed &&\
bash src/date.sh results/coverage/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_sort_in_genes.bed 
