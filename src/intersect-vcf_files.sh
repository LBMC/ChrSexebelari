# Intersect 2 raw VCF files from male and female variants calling 
# Make summary of intersected and variantspresent in one sexe from both files 

# done on local computer 

# If VCF are not VCF.gz: use "bgzip -c input.vcf > output.vcf.gz"

# Index each VCF
bcftools index results/call_var/2017_09_21_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz &&\
bcftools index results/call_var/2017_09_21_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz &&\

# Make intersection: this generates 2 vcf files: 
# - 0000.vcf which represents variants in MRDR6 (1st sample) and not in MRDR5 (2nd sample)
# - 0002.vcf which represents variants present in both pools
bcftools isec results/call_var/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz results/call_var/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz -p results/call_var/inter_female_male_realign_indels -w 1 &&\

# Make summary with rtg
bin/rtg vcfstats results/call_var/inter_female_male_realign_indels/0000.vcf > results/call_var/inter_female_male_realign_indels/summary_in_MRDR6_not_in_MRDR5_0000.txt &&\
bin/rtg vcfstats results/call_var/inter_female_male_realign_indels/0002.vcf > results/call_var/inter_female_male_realign_indels/summary_in_both_0002.txt 
