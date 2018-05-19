bwa index 2017_09_08_combined_belari_coli.fasta

bwa mem -x ont2d data/reference_genomes/2017_09_08_combined_belari_coli.fasta data/nanopore_merged_files/2018-05-17-Nanopore_all_reads_best_qual.fastq > results/nanopore_against_ecoli_genome/2018-05-18-Nanopore_ecoli_align.bam

samtools view -f4 results/nanopore_against_ecoli_genome/2018-05-18-Nanopore_ecoli_align.bam > results/nanopore_against_ecoli_genome/2018-05-19-Nanopore_ecoli_align.unmapped.bam

cut -f1 results/nanopore_against_ecoli_genome/2018-05-19-Nanopore_ecoli_align.unmapped.bam | sort | uniq > results/nanopore_against_ecoli_genome/2018-05-19-Nanopore_ids_unmapped_vs_ecoli.lst

seqtk subseq data/nanopore_merged_files/2018-05-17-Nanopore_all_reads_best_qual.fastq results/nanopore_against_ecoli_genome/2018-05-19-Nanopore_ids_unmapped_vs_ecoli.lst > results/nanopore_against_ecoli_genome/2018-05-17-Nanopore_all_reads_best_qual_wo_ecoli.fastq

