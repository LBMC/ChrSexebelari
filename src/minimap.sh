query=
subject=
output=

~/Programs/minimap2/minimap2 -ax map-ont ~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta ~/Documents/stage_mbelari/results/nanopore_against_ecoli_genome/2018-05-24-Nanopore_all_reads_best_qual_wo_ecoli.fastq > 2018-06-18-Nanopore_reads_vs_hybrid.sam

samtools view -Sb 2018-06-18-Nanopore_reads_vs_hybrid.sam > 2018-06-18-Nanopore_reads_vs_hybrid.bam

samtools sort 2018-06-18-Nanopore_reads_vs_hybrid.bam -o 2018-06-18-Nanopore_reads_vs_hybrid.sorted.bam

samtools index 2018-06-18-Nanopore_reads_vs_hybrid.sorted.bam


