bwa mem -x ont2d ~/Documents/stage_mbelari/data/reference_genomes/ADBT01.1.fsa_nt.gz ~/Documents/stage_mbelari/data/nanopore_merged_files/2018-05-24-Nanopore_all_reads_best_qual.fastq > ~/Documents/stage_mbelari/results/nanopore_against_ecoli_genome/2018-05-24-Nanopore_ecoli_align.bam

samtools view -f4 ~/Documents/stage_mbelari/results/nanopore_against_ecoli_genome/2018-05-24-Nanopore_ecoli_align.bam > ~/Documents/stage_mbelari/results/nanopore_against_ecoli_genome/2018-05-24-Nanopore_ecoli_align.unmapped.bam

samtools view -F4 ~/Documents/stage_mbelari/results/nanopore_against_ecoli_genome/2018-05-24-Nanopore_ecoli_align.bam > ~/Documents/stage_mbelari/results/nanopore_against_ecoli_genome/2018-05-24-Nanopore_ecoli_align.mapped.bam

cut -f1 ~/Documents/stage_mbelari/results/nanopore_against_ecoli_genome/2018-05-24-Nanopore_ecoli_align.unmapped.bam | sort | uniq > ~/Documents/stage_mbelari/results/nanopore_against_ecoli_genome/2018-05-24-Nanopore_ids_unmapped_vs_ecoli.lst

seqtk subseq ~/Documents/stage_mbelari/data/nanopore_merged_files/2018-05-24-Nanopore_all_reads_best_qual.fastq ~/Documents/stage_mbelari/results/nanopore_against_ecoli_genome/2018-05-24-Nanopore_ids_unmapped_vs_ecoli.lst > ~/Documents/stage_mbelari/results/nanopore_against_ecoli_genome/2018-05-24-Nanopore_all_reads_best_qual_wo_ecoli.fastq

