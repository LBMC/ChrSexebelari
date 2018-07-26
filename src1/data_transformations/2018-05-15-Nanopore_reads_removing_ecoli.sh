bwa index data/reference_genomes/ADBT01.1.fsa_nt.gz

bwa mem -x ont2d data/reference_genomes/ADBT01.1.fsa_nt.gz data/nanopore_merged_files/2018-05-17-Nanopore_all_reads_best_qual.fastq > results/nanopore_against_ecoli_genome/2018-05-18-Nanopore_ecoli_align.bam

samtools view -f4 results/nanopore_against_ecoli_genome/2018-05-18-Nanopore_ecoli_align.bam > results/nanopore_against_ecoli_genome/2018-05-19-Nanopore_ecoli_align.unmapped.bam

samtools view -F4 results/nanopore_against_ecoli_genome/2018-05-18-Nanopore_ecoli_align.bam > results/nanopore_against_ecoli_genome/2018-05-19-Nanopore_ecoli_align.mapped.bam

cut -f1 results/nanopore_against_ecoli_genome/2018-05-19-Nanopore_ecoli_align.unmapped.bam | sort | uniq > results/nanopore_against_ecoli_genome/2018-05-19-Nanopore_ids_unmapped_vs_ecoli.lst

seqtk subseq data/nanopore_merged_files/2018-05-17-Nanopore_all_reads_best_qual.fastq results/nanopore_against_ecoli_genome/2018-05-19-Nanopore_ids_unmapped_vs_ecoli.lst > results/nanopore_against_ecoli_genome/2018-05-19-Nanopore_all_reads_best_qual_wo_ecoli.fastq

