# Fastq to fasta

cat results/mapping/mapped/2017_10_30_MRDR6_trim_Mbelari_Ecoli_unmapped.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > results/mapping/mapped/2017_11_06_MRDR6_trim_Mbelari_Ecoli_unmapped.fa &&\

cat results/mapping/mapped/2017_10_30_MRDR5_trim_Mbelari_Ecoli_unmapped.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > results/mapping/mapped/2017_11_06_MRDR5_trim_Mbelari_Ecoli_unmapped.fa
