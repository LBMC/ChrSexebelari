assembly=~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta
fasta_file=~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-24-Genes_absent_in_females.fasta

outputfile=2018-07-24-BLAST_Genes_absent_in_females.txt
outputfilealign=2018-07-24-BLAST_Genes_absent_in_females.align.txt

#makeblastdb -in ${assembly} -parse_seqids -dbtype nucl

blastn -query ${fasta_file} -db ${assembly} -out ${outputfile} -outfmt 6
blastn -query ${fasta_file} -db ${assembly} -out ${outputfilealign} -outfmt 4
