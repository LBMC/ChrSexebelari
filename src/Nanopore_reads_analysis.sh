INPUTreads=~/Documents/stage_mbelari/results/nanopore_against_ecoli_genome/2018-05-24-Nanopore_all_reads_best_qual_wo_ecoli.fastq
OUTdir=~/Documents/stage_mbelari/results/nanopore_stats

NanoStat --fastq ${INPUTreads} --outdir ${OUTdir} -p 2018-07-07-Nanopore_clean -n Nanopore_clean
