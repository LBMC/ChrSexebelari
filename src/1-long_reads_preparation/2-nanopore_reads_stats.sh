fastq=results/cleaned_from_ecoli/2018-05-24-Nanopore_all_reads_best_qual_wo_ecoli.fastq

NanoStat --fastq ${fastq} -o results/nanopore_stats -n Nanopore_clean
