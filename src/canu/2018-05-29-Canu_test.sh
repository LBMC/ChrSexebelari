export PATH=$PATH:~/Programs/canu-1.7/Linux-amd64/bin

canu -d Nanopore_reads_new_subset -p Nanopore_reads_new_subset genomeSize=57k minReadLength=500 minOverlapLength=250 -nanopore-raw ~/Documents/stage_mbelari/results/nanopore_subset/2018-05-29-Nanopore_MBELA01182_subset_2.fq

