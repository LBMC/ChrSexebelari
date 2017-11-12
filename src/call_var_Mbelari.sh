#samtools 1.3
#bcftools 1.3

# To call variants: samtools mpileup -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -ABQ0 -L 100000 -uf $1 $2 | bcftools call -mv --output-type z -p 0.1 -o $3
samtools mpileup -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -ABQ0 -L 100000 -uf data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_scaffolds_repeatmasker_masked.fa results/mapping/without_duplicates/MRDR5_trim_Mbelari_mapped_rmdup.bam | bcftools call -mv --output-type z -p 0.1 -o results/call_var/MRDR5_trim_Mbelari_mapped_rmdup.vcf.gz &&\

samtools mpileup -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -ABQ0 -L 100000 -uf data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_scaffolds_repeatmasker_masked.fa results/mapping/without_duplicates/MRDR6_trim_Mbelari_mapped_rmdup.bam | bcftools call -mv --output-type z -p 0.1 -o results/call_var/MRDR6_trim_Mbelari_mapped_rmdup.vcf.gz &&\

# To call variants (-mv) and get pileup in VCF format (-mA) after realignment around INDELs: samtools mpileup -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -ABQ0 -L 100000 -uf $1 $2 | bcftools call -mA --output-type z -p 0.1 -o $3
samtools mpileup -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -ABQ0 -L 100000 -uf data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_scaffolds_repeatmasker_masked.fa results/mapping/without_duplicates/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.bam | bcftools call -mv --output-type z -p 0.1 -o results/call_var/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz &&\

samtools mpileup -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -ABQ0 -L 100000 -uf data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_scaffolds_repeatmasker_masked.fa results/mapping/without_duplicates/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.bam | bcftools call -mv --output-type z -p 0.1 -o results/call_var/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz &&\


samtools mpileup -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -ABQ0 -L 100000 -uf data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_scaffolds_repeatmasker_masked.fa results/mapping/without_duplicates/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.bam | bcftools call -mA --output-type z -p 0.1 -o results/call_var/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_A.vcf.gz &&\

samtools mpileup -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -ABQ0 -L 100000 -uf data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_scaffolds_repeatmasker_masked.fa results/mapping/without_duplicates/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.bam | bcftools call -mA --output-type z -p 0.1 -o results/call_var/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_A.vcf.gz 

