./bin/vcftools --gzvcf results/call_var/MRDR5_trim_Mbelari_mapped_rmdup.vcf.gz --keep-only-indels --recode --recode-INFO-all | gzip -c > results/call_var/MRDR5_trim_Mbelari_mapped_rmdup_indels.vcf.gz &&\

./bin/vcftools --gzvcf results/call_var/MRDR6_trim_Mbelari_mapped_rmdup.vcf.gz --keep-only-indels --recode --recode-INFO-all | gzip -c > results/call_var/MRDR6_trim_Mbelari_mapped_rmdup_indels.vcf.gz 

#./bin/vcftools --gvcf results/call_var/2017_09_13_MRDR6_trim_Mbelari_mapped_sort_rmdup.vcf --keep-only-indels --recode-INFO-all --out results/call_var/MRDR6_trim_Mbelari_mapped_sort_rmdup_indels --recode 
