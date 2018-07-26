#Convert Fastq a fasta
FastqINPUT=results/cleaned_from_ecoli/2018-05-24-Nanopore_all_reads_best_qual_wo_ecoli.fastq
FastaOUTPUT=results/cleaned_from_ecoli/2018-05-24-Nanopore_all_reads_best_qual_wo_ecoli.fasta
seqtk seq -a ${FastqINPUT} > ${FastaOUTPUT}
