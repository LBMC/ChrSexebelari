# steps of the analysis of the repreat content

# female_JU2817
mkdir -p results/dnaPipeTE/female_JU2817
python3 src/dnaPipeTE/dnaPipeTE.py \
-input results/mapping/mapped/2017_09_13_MRDR5_trim_Mbelari_mapped_sort.fastq.gz \
-output results/dnaPipeTE/female_JU2817 \
-cpu 10 \
-genome_size 162071050 \
-genome_coverage 43.7 \
-sample_number 3 2> results/female_JU2817_report.txt

# male_JU2817
mkdir -p results/dnaPipeTE/male_JU2817
python3 src/dnaPipeTE/dnaPipeTE.py \
-input results/mapping/mapped/2017_09_13_MRDR6_trim_Mbelari_mapped_sort.fastq.gz \
-output results/dnaPipeTE/male_JU2817 \
-cpu 10 \
-genome_size 162071050 \
-genome_coverage 24.8 \
-sample_number 3 2> results/male_JU2817_report.txt
