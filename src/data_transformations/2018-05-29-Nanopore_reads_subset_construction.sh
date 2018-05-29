#bwa index ~/Documents/stage_mbelari/data/sample/MBELA01182_RNAP2.fa 

#bwa mem -x ont2d ~/Documents/stage_mbelari/data/sample/MBELA01182_RNAP2.fa ~/Documents/stage_mbelari/data/sample/MBELA01182subset.fastq > ~/Documents/stage_mbelari/results/nanopore_subset/2018-05-29-Nanopore_MBELA01182_subset.bam

#samtools sort ~/Documents/stage_mbelari/results/nanopore_subset/2018-05-29-Nanopore_MBELA01182_subset.bam -o ~/Documents/stage_mbelari/results/nanopore_subset/2018-05-29-Nanopore_MBELA01182_subset.sorted.bam

#samtools index ~/Documents/stage_mbelari/results/nanopore_subset/2018-05-29-Nanopore_MBELA01182_subset.sorted.bam

#samtools faidx ~/Documents/stage_mbelari/data/sample/MBELA01182_RNAP2.fa



bwa index ~/Documents/stage_mbelari/data/sample/MBELA01182_RNAP2.fa 

bwa mem -x ont2d ~/Documents/stage_mbelari/data/sample/MBELA01182_RNAP2.fa ~/Documents/stage_mbelari/results/nanopore_against_ecoli_genome/2018-05-24-Nanopore_all_reads_best_qual_wo_ecoli.fastq > ~/Documents/stage_mbelari/results/nanopore_subset/2018-05-29-Nanopore_MBELA01182_subset_2.bam

samtools view -F4 ~/Documents/stage_mbelari/results/nanopore_subset/2018-05-29-Nanopore_MBELA01182_subset_2.bam > ~/Documents/stage_mbelari/results/nanopore_subset/2018-05-29-Nanopore_MBELA01182_subset_2.mapped.bam

cut -f1 ~/Documents/stage_mbelari/results/nanopore_subset/2018-05-29-Nanopore_MBELA01182_subset_2.mapped.bam | sort | uniq > ~/Documents/stage_mbelari/results/nanopore_subset/2018-05-29-Nanopore_MBELA01182_subset_2_mapped.lst
