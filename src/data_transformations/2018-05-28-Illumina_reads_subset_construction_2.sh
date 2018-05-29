bowtie2 -x ~/Documents/stage_mbelari/results/illumina_subset/MBELA01182-bowtie2 -1 ~/Documents/stage_mbelari/data/sample/MR_350_mapped_to_Mbelari-01182_1.fastq -2 ~/Documents/stage_mbelari/data/sample/MR_350_mapped_to_Mbelari-01182_2.fastq -S ~/Documents/stage_mbelari/results/illumina_subset/2018-05-28-Illumina_subset_MBELA01182.sam

samtools view -Sb ~/Documents/stage_mbelari/results/illumina_subset/2018-05-28-Illumina_subset_MBELA01182.sam > ~/Documents/stage_mbelari/results/illumina_subset/2018-05-28-Illumina_subset_MBELA01182.bam

samtools sort ~/Documents/stage_mbelari/results/illumina_subset/2018-05-28-Illumina_subset_MBELA01182.bam -o ~/Documents/stage_mbelari/results/illumina_subset/2018-05-28-Illumina_subset_MBELA01182.sorted.bam

samtools index ~/Documents/stage_mbelari/results/illumina_subset/2018-05-28-Illumina_subset_MBELA01182.sorted.bam

samtools faidx ~/Documents/stage_mbelari/data/sample/MBELA01182_RNAP2.fa
