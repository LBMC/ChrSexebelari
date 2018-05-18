bwa index 2017_09_08_combined_belari_coli.fasta

bwa mem -x ont2d ~/Documents/stage_mbelari/data/reference_genomes/2017_09_08_combined_belari_coli.fasta ~/Documents/stage_mbelari/data/nanopore_merged_files/2018-05-17-Nanopore_all_reads_best_qual.fastq > ~/Documents/stage_mbelari/data/nanopore_merged_files/2018-05-18-Nanopore_ecoli_align.bam


