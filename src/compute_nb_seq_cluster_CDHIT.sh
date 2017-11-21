TotalClusters="$(grep -c '>Cluster' results/mapping/mapped/2017_11_21_cdhit.fa.clstr)"

awk -F" " 'BEGIN {print "cluster" "\t" "nb"} $1 ~ /^>/ {if (count) print $2 "\t" count-1; count=0; next} {count++} END {print "'$TotalClusters'" "\t" count-1}' results/mapping/mapped/2017_11_21_cdhit.fa.clstr > results/mapping/mapped/2017_11_22_cdhit.fa.clstr.size

