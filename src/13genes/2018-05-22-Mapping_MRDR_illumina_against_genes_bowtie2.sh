bowtie2-build data/raw/2018-05-18-Complete_13_genes.fasta results/13genes/13_genes

bowtie2 -x results/13genes/13_genes -U data/raw/2017_09_13_MRDR5_trim_Mbelari_mapped_sort.fastq.gz -S results/13genes/2018-05-22-MRDR5_Mbelai_vs_13_genes.sam --no-unal

bowtie2 -x results/13genes/13_genes -U data/raw/2017_09_13_MRDR6_trim_Mbelari_mapped_sort.fastq.gz -S results/13genes/2018-05-22-MRDR6_Mbelai_vs_13_genes.sam --no-unal

samtools view -Sb results/13genes/2018-05-22-MRDR5_Mbelai_vs_13_genes.sam > results/13genes/2018-05-22-MRDR5_Mbelai_vs_13_genes.bam

samtools view -Sb results/13genes/2018-05-22-MRDR6_Mbelai_vs_13_genes.sam > results/13genes/2018-05-22-MRDR6_Mbelai_vs_13_genes.bam

samtools sort results/13genes/2018-05-22-MRDR5_Mbelai_vs_13_genes.bam -o results/13genes/2018-05-22-MRDR5_Mbelai_vs_13_genes_sorted.bam

samtools sort results/13genes/2018-05-22-MRDR6_Mbelai_vs_13_genes.bam -o results/13genes/2018-05-22-MRDR6_Mbelai_vs_13_genes_sorted.bam

samtools index results/13genes/2018-05-22-MRDR5_Mbelai_vs_13_genes_sorted.bam

samtools index results/13genes/2018-05-22-MRDR6_Mbelai_vs_13_genes_sorted.bam

samtools faidx data/raw/2018-05-18-Complete_13_genes.fasta
