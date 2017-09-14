# steps of the analysis of the repreat content

# female_JU2817
mkdir -p results/dnaPipeTE/female_JU2817
python3 src/dnaPipeTE/dnaPipeTE.py \
-input data/fastq/2017_09_13_female_JU2817_R1.fastq.gz \
-output results/dnaPipeTE/female_JU2817 \
-cpu 12 \
-genome_size 162071050 \
-genome_coverage 15 \
-sample_number 5

# male_JU2817
mkdir -p results/dnaPipeTE/male_JU2817
python3 src/dnaPipeTE/dnaPipeTE.py \
-input data/fastq/2017_09_13_male_JU2817_R1.fastq.gz \
-output results/dnaPipeTE/male_JU2817 \
-cpu 12 \
-genome_size 162071050 \
-genome_coverage 15 \
-sample_number 5
