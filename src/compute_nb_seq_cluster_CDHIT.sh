# Parser to extract cluster sizes from CDHIT output

TotalClusters="$(grep -c '>Cluster' results/mapping/unmapped/2017_11_21_cdhit.fa.clstr)"

awk -F" " 'BEGIN {print "cluster" "\t" "nb"} $1 ~ /^>/ {if (count) print $2 "\t" count-1; count=0; next} {count++} END {print "'$TotalClusters'" "\t" count-1}' results/mapping/unmapped/2017_11_21_cdhit.fa.clstr > results/mapping/unmapped/2017_11_22_cdhit.fa.clstr.size

