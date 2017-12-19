# Create intervals to realign around INDELs with GATK 
# Note that the input bam must have been indexed
# done on local computer

# Create intervals
java -jar bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R data/ReferenceGenomes/2017_09_08_combined_belari_coli.fasta -I results/mapping/without_duplicates/2017_09_20_MRDR5_trim_Mbelari_mapped_rmdup_rg.bam -o results/mapping/without_duplicates/2017_09_20_MRDR5_trim_Mbelari_mapped_rmdup_rg.intervals &&\

java -jar bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R data/ReferenceGenomes/2017_09_08_combined_belari_coli.fasta -I results/mapping/without_duplicates/2017_09_20_MRDR6_trim_Mbelari_mapped_rmdup_rg.bam -o results/mapping/without_duplicates/2017_09_20_MRDR6_trim_Mbelari_mapped_rmdup_rg.intervals &&\

# Do the realignment
java -jar bin/GenomeAnalysisTK.jar -T IndelRealigner -I results/mapping/without_duplicates/2017_09_20_MRDR5_trim_Mbelari_mapped_rmdup_rg.bam -R data/ReferenceGenomes/2017_09_08_combined_belari_coli.fasta -targetIntervals results/mapping/without_duplicates/2017_09_20_MRDR5_trim_Mbelari_mapped_rmdup_rg.intervals -o results/mapping/without_duplicates/2017_09_20_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.bam &&\

java -jar bin/GenomeAnalysisTK.jar -T IndelRealigner -I results/mapping/without_duplicates/2017_09_20_MRDR6_trim_Mbelari_mapped_rmdup_rg.bam -R data/ReferenceGenomes/2017_09_08_combined_belari_coli.fasta -targetIntervals results/mapping/without_duplicates/2017_09_20_RDR6_trim_Mbelari_mapped_rmdup_rg.intervals -o results/mapping/without_duplicates/2017_09_20_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.bam &&\

# Index output
samtools index results/mapping/without_duplicates/2017_09_20_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.bam &&\
samtools index results/mapping/without_duplicates/2017_09_20_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.bam

