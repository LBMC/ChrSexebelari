awk 'BEGIN{-F""} 
$0 ~ /^>/ {print $0} 
END {print "name of seq"}' $1 > $2
