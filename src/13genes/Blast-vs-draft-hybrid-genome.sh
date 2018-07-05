assembly=~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta
fasta_file=~/Documents/stage_mbelari/results/13genes/version2/sex-determination-genes.fasta

outputfile=2018-07-04-sex-deter-genes.txt
outputfilealign=2018-07-sex-deter-genes.txt

#makeblastdb -in ${assembly} -parse_seqids -dbtype nucl

blastn -query ${fasta_file} -db ${assembly} -out ${outputfile} -outfmt 6
blastn -query ${fasta_file} -db ${assembly} -out ${outputfilealign}
