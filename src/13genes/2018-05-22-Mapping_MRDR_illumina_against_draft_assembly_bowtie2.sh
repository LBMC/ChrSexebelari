bowtie2-build data/raw/mbela_assembly_JU2817.fa results/13genes/2018-05-22-mbela_assembly_JU2817

bowtie2 -x results/13genes/2018-05-22-mbela_assembly_JU2817 -U data/raw/2017_09_13_MRDR5_trim_Mbelari_mapped_sort.fastq.gz -S results/13genes/2018-05-22-MRDR5_vs_Mbelari_assembly.sam --no-unal

bowtie2 -x results/13genes/2018-05-22-mbela_assembly_JU2817 -U data/raw/2017_09_13_MRDR6_trim_Mbelari_mapped_sort.fastq.gz -S results/13genes/2018-05-22-MRDR6_vs_Mbelari_assembly.sam --no-unal

samtools view -Sb results/13genes/2018-05-22-MRDR5_vs_Mbelari_assembly.sam > results/13genes/2018-05-22-MRDR5_vs_Mbelari_assembly.bam

samtools view -Sb results/13genes/2018-05-22-MRDR6_vs_Mbelari_assembly.sam > results/13genes/2018-05-22-MRDR6_vs_Mbelari_assembly.bam

samtools sort results/13genes/2018-05-22-MRDR5_vs_Mbelari_assembly.bam -o results/13genes/2018-05-22-MRDR5_vs_Mbelari_assembly.sorted.bam

samtools sort results/13genes/2018-05-22-MRDR6_vs_Mbelari_assembly.bam -o results/13genes/2018-05-22-MRDR6_vs_Mbelari_assembly.sorted.bam

samtools index results/13genes/2018-05-22-MRDR5_vs_Mbelari_assembly.sorted.bam

samtools index results/13genes/2018-05-22-MRDR6_vs_Mbelari_assembly.sorted.bam

samtools faidx data/raw/mbela_assembly_JU2817.fa
