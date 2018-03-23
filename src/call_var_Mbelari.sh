# Call variants before and after realignment around INDELs
# Done locally

# before realignment around INDELS
samtools mpileup -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -ABQ0 -L 100000 -uf data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_scaffolds_repeatmasker_masked.fa results/mapping/without_duplicates/2017_09_20_MRDR5_trim_Mbelari_mapped_rmdup.bam | bcftools call -mv --output-type z -p 0.1 -o results/call_var/2017_09_20_MRDR5_trim_Mbelari_mapped_rmdup.vcf.gz &&\

samtools mpileup -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -ABQ0 -L 100000 -uf data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_scaffolds_repeatmasker_masked.fa results/mapping/without_duplicates/2017_09_20_MRDR6_trim_Mbelari_mapped_rmdup.bam | bcftools call -mv --output-type z -p 0.1 -o results/call_var/2017_09_20_MRDR6_trim_Mbelari_mapped_rmdup.vcf.gz &&\

# after realignment around INDELS
samtools mpileup -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -ABQ0 -L 100000 -uf data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_scaffolds_repeatmasker_masked.fa results/mapping/without_duplicates/2017_09_20_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.bam | bcftools call -mv --output-type z -p 0.1 -o results/call_var/2017_09_21_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz &&\

samtools mpileup -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -ABQ0 -L 100000 -uf data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_scaffolds_repeatmasker_masked.fa results/mapping/without_duplicates/2017_09_20_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.bam | bcftools call -mv --output-type z -p 0.1 -o results/call_var/2017_09_21_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz 

