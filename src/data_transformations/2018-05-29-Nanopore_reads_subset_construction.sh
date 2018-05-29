bwa index ~/Documents/stage_mbelari/data/sample/MBELA01182_RNAP2.fa 

bwa mem -x ont2d ~/Documents/stage_mbelari/data/sample/MBELA01182_RNAP2.fa ~/Documents/stage_mbelari/data/sample/MBELA01182subset.fastq > ~/Documents/stage_mbelari/results/nanopore_against_ecoli_genome/2018-05-18-Nanopore_MBELA01182_subset.bam

samtools sort ~/Documents/stage_mbelari/results/illumina_subset/2018-05-28-Illumina_subset_MBELA01182.bam -o ~/Documents/stage_mbelari/results/illumina_subset/2018-05-28-Illumina_subset_MBELA01182.sorted.bam

samtools index ~/Documents/stage_mbelari/results/illumina_subset/2018-05-28-Illumina_subset_MBELA01182.sorted.bam

samtools faidx ~/Documents/stage_mbelari/data/sample/MBELA01182_RNAP2.fa
