working_directory=$1

cd ${working_directory}

bowtie2-build /scratch/lestrada/stage_mbelari/results/hybrid_test/full_test/test1/consensus_dir/final_assembly.fasta hybrid_index

bowtie2 -x hybrid_index -U /scratch/lestrada/stage_mbelari/data/raw/2017_09_13_MRDR5_trim_Mbelari_mapped_sort.fastq.gz -S 2018-06-05-MRDR5_vs_Hybrid_assembly.sam

bowtie2 -x hybrid_index -U /scratch/lestrada/stage_mbelari/data/raw/2017_09_13_MRDR6_trim_Mbelari_mapped_sort.fastq.gz -S 2018-06-05-MRDR6_vs_Hybrid_assembly.sam

samtools view -Sb 2018-06-05-MRDR5_vs_Hybrid_assembly.sam > 2018-06-05-MRDR5_vs_Hybrid_assembly.bam

samtools view -Sb 2018-06-05-MRDR6_vs_Hybrid_assembly.sam > 2018-06-05-MRDR6_vs_Hybrid_assembly.bam

samtools sort 2018-06-05-MRDR5_vs_Hybrid_assembly.bam -o 2018-06-05-MRDR5_vs_Hybrid_assembly.sorted.bam

samtools sort 2018-06-05-MRDR6_vs_Hybrid_assembly.bam -o 2018-06-05-MRDR6_vs_Hybrid_assembly.sorted.bam

samtools index 2018-06-05-MRDR5_vs_Hybrid_assembly.sorted.bam

samtools index 2018-06-05-MRDR6_vs_Hybrid_assembly.sorted.bam

#samtools faidx ~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta

