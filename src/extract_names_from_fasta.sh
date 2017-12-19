# Extract names of sequences from a fasta.
# done on local computer 
# $1 is input file and $2 is output for saving sequence name from fasta

awk 'BEGIN{-F""} 
$0 ~ /^>/ {print $0} 
END {print "name of seq"}' $1 > $2
