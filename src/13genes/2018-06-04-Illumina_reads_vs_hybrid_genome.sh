bowtie2-build ~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta hybrid_index

bowtie2 --very-sensitive -x ~/Documents/stage_mbelari/results/13genes/hybrid_index -1 ~/Documents/stage_mbelari/data/raw/MR_350_clean_1.fastq -2 ~/Documents/stage_mbelari/data/raw/MR_350_clean_2.fastq -S Illumina_reads_vs_hgenome

samtools view -Sb ~/Documents/stage_mbelari/results/13genes/Illumina_reads_vs_hgenome.sam > ~/Documents/stage_mbelari/results/13genes/Illumina_reads_vs_hgenome.bam

samtools view -f 2 ~/Documents/stage_mbelari/results/13genes/Illumina_reads_vs_hgenome.bam > ~/Documents/stage_mbelari/results/13genes/Illumina_reads_vs_hgenome.paired.bam

samtools sort ~/Documents/stage_mbelari/results/13genes/Illumina_reads_vs_hgenome.paired.bam > ~/Documents/stage_mbelari/results/13genes/Illumina_reads_vs_hgenome.sorted.bam