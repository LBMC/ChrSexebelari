# Summarize the number of variants found in vcf files from RTG
# done locally

##### raw variants 
bin/rtg vcfstats results/call_var/2017_09_20_MRDR5_trim_Mbelari_mapped_rmdup.vcf.gz > results/call_var/summary_raw_var_female_JU2917.txt &&\

bin/rtg vcfstats results/call_var/2017_09_20_MRDR6_trim_Mbelari_mapped_rmdup.vcf.gz > results/call_var/summary_raw_var_male_JU2917.txt &&\

##### after realignment around indels
bin/rtg vcfstats results/call_var/2017_09_21_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz > results/call_var/summary_realign_indels_var_female_JU2917.txt &&\

bin/rtg vcfstats results/call_var/2017_09_21_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz > results/call_var/summary_realign_indels_var_male_JU2917.txt


