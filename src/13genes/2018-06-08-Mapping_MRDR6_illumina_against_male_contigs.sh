
bowtie2-build ~/Documents/stage_mbelari/results/13genes/2018-06-08-Contigs_males.fasta male_contigs

bowtie2 -x male_contigs -U ~/Documents/stage_mbelari/data/raw/2017_09_13_MRDR6_trim_Mbelari_mapped_sort.fastq.gz -S 2018-06-08-MRDR6_vs_contigs_males.sam

samtools view -Sb 2018-06-08-MRDR6_vs_contigs_males.sam > 2018-06-08-MRDR6_vs_contigs_males.bam

samtools sort 2018-06-08-MRDR6_vs_contigs_males.bam -o 2018-06-08-MRDR6_vs_contigs_males.sorted.bam

samtools index 2018-06-08-MRDR6_vs_contigs_males.sorted.bam

#samtools faidx ~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta

