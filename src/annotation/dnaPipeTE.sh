INPUT=~/stage_mbelari/data/2017_09_14_MRDR5_trim_Mbelari_mapped_sort_merged.fastq
OUTPUT=~/stage_mbelari/results/dnaPipeTE/mbelari_MRDR5_01

python3 /opt/dnaPipeTE/dnaPipeTE.py -input ${INPUT} -output ${OUTPUT} -genome_size 129550000 -genome_coverage 0.1 -keep_Trinity_output 
