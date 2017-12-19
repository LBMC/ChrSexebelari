# Extract first element of names of sequences from a fasta.
# done on local computer 

awk 'BEGIN{-F" "} $0 ~ /^>/ {print substr($1,2)} END {}' 2017_09_12_contigs_name_Ecoli.txt > 2017_09_12_contigs_short_name_Ecoli.txt &&\

awk 'BEGIN{-F" "} $0 ~ /^>/ {print substr($1,2)} END {}' 2017_09_12_contigs_name_Mbelari.txt > 2017_09_12_contigs_short_name_Mbelari.txt


