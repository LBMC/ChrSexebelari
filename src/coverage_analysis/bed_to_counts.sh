#The bam files were filtered by quality

INPUTdeb=~/Documents/stage_mbelari/results/coverage_analysis/2018_06_21_MRDR5_vs_Hybrid_assembly.sorted.filtered.bed
GENOMEsizes=~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt
OUTPUT=~/Documents/stage_mbelari/results/coverage_analysis/2018_06_21_Illumina_vs_Hybrid_assembly.sorted.filtered.bed

bedtools genomecov -i $INPUTdeb -g $GENOMEsizes > $OUTPUT

