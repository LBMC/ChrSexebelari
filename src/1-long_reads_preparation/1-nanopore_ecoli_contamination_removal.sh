cat data/raw/FC17_all_reads_best_qual.fastq data/raw/FC18_all_reads_best_qual.fastq data/raw/FC19_all_reads_best_qual.fastq data/raw/FC20_all_reads_best_qual.fastq > results/cleaned_from_ecoli/2018-05-24-Nanopore_all_reads_best_qual.fastq

bwa index data/reference_genomes/ADBT01.1.fsa_nt.gz

bwa mem -x ont2d data/reference_genomes/ADBT01.1.fsa_nt.gz results/cleaned_from_ecoli/2018-05-24-Nanopore_all_reads_best_qual.fastq > results/cleaned_from_ecoli/2018-05-24-Nanopore_ecoli_align.bam

samtools view -f4 results/cleaned_from_ecoli/2018-05-24-Nanopore_ecoli_align.bam > results/cleaned_from_ecoli/2018-05-24-Nanopore_ecoli_align.unmapped.bam

samtools view -F4 results/cleaned_from_ecoli/2018-05-24-Nanopore_ecoli_align.bam > results/cleaned_from_ecoli/2018-05-24-Nanopore_ecoli_align.mapped.bam

cut -f1 results/cleaned_from_ecoli/2018-05-24-Nanopore_ecoli_align.unmapped.bam | sort | uniq > results/cleaned_from_ecoli/2018-05-24-Nanopore_ids_unmapped_vs_ecoli.lst

seqtk subseq results/cleaned_from_ecoli/2018-05-24-Nanopore_all_reads_best_qual.fastq results/cleaned_from_ecoli/2018-05-24-Nanopore_ids_unmapped_vs_ecoli.lst > results/cleaned_from_ecoli/2018-05-24-Nanopore_all_reads_best_qual_wo_ecoli.fastq

