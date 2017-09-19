# steps of the analysis of the repreat content

# female_JU2817
mkdir -p results/dnaPipeTE/female_JU2817
python3 src/dnaPipeTE/dnaPipeTE.py \
-input results/mapping/mapped/2017_09_13_MRDR5_trim_Mbelari_mapped_sort.fastq.gz \
-output results/dnaPipeTE/female_JU2817 \
-cpu 12 \
-genome_size 162071050 \
-genome_coverage 22.9 \
-sample_number 3

# male_JU2817
mkdir -p results/dnaPipeTE/male_JU2817
python3 src/dnaPipeTE/dnaPipeTE.py \
-input results/mapping/mapped/2017_09_13_male_JU2817_R1.fastq.gz \
-output results/dnaPipeTE/male_JU2817 \
-cpu 12 \
-genome_size 162071050 \
-genome_coverage 13.6 \
-sample_number 5
