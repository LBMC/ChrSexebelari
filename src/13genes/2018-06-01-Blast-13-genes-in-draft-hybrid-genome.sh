assembly=~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta
fasta_file=~/Documents/stage_mbelari/results/13genes/version2/2018-06-20-Top_50_genes_fem_abs.fasta

outputfile=2018-06-20-Genes_absent_females.txt
outputfilealign=2018-06-20-Genes_absent_females-align.txt

#makeblastdb -in ${assembly} -parse_seqids -dbtype nucl

blastn -query ${fasta_file} -db ${assembly} -out ${outputfile} -outfmt 6
blastn -query ${fasta_file} -db ${assembly} -out ${outputfilealign}
