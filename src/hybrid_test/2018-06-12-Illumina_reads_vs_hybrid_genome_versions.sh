assembly= ~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta
pe1= ~/Documents/stage_mbelari/data/raw/MR_350_clean_1.fastq
pe2= ~/Documents/stage_mbelari/data/raw/MR_350_clean_2.fastq

bowtie2-build ${assembly} hybrid_index

bowtie2 --very-sensitive -x hybrid_index -1 ${pe1} -2 ${pe2} -S 2018-06-14-Illumina_reads_vs_hgenome_3_250.sam

samtools view -Sb 2018-06-14-Illumina_reads_vs_hgenome_3_250.sam > 2018-06-14-Illumina_reads_vs_hgenome_3_250.bam

samtools sort 2018-06-14-Illumina_reads_vs_hgenome_3_250.bam > 2018-06-14-Illumina_reads_vs_hgenome_3_250.sorted.bam

samtools index 2018-06-14-Illumina_reads_vs_hgenome_3_250.sorted.bam

#samtools faidx ~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta

