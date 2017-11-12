# Intersect 2 vcf files 

# if not done, use bgzip -c input.vcf > output.vcf.gz (

#bcftools 1.3 and bgzip 1.1

bcftools index results/call_var/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz &&\
bcftools index results/call_var/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz &&\

bcftools isec results/call_var/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz results/call_var/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz -p results/call_var/inter_female_male_realign_indels -w 1 &&\

bin/rtg vcfstats results/call_var/inter_female_male_realign_indels/0000.vcf > results/call_var/inter_female_male_realign_indels/summary_in_MRDR6_not_in_MRDR5_0000.txt &&\
bin/rtg vcfstats results/call_var/inter_female_male_realign_indels/0002.vcf > results/call_var/inter_female_male_realign_indels/summary_in_both_0002.txt 
