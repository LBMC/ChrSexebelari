awk 'BEGIN{-F""} 
$0 ~ /^>/ {print $0} 
END {print "name of contigs"}' $1 > $2
