#The bam files were filtered by quality

INPUTcov=results/coverage_analysis/2018-06-12-MRDR6_vs_Hybrid_assembly.sorted.filtered.bam
OUTPUTcov=results/coverage_analysis/2018_07_25_MRDR6_vs_Hybrid_assembly.sorted.filtered.bed

bedtools genomecov -d -ibam $INPUTcov > $OUTPUTcov

INPUTcov=results/coverage_analysis/2018-06-12-MRDR5_vs_Hybrid_assembly.sorted.filtered.bam
OUTPUTcov=results/coverage_analysis/2018_07_25_MRDR5_vs_Hybrid_assembly.sorted.filtered.bed

bedtools genomecov -d -ibam $INPUTcov > $OUTPUTcov


BED1=results/coverage_analysis/2018_07_25_MRDR6_vs_Hybrid_assembly.sorted.filtered.bed
BED2=results/coverage_analysis/2018_07_25_MRDR5_vs_Hybrid_assembly.sorted.filtered.bed
OUTPUTcov=results/coverage_analysis/2018_07_25_MRDR5_MRDR6_vs_Hybrid_assembly.sorted.filtered.bed

paste $BED1 $BED2 | awk 'BEGIN{print "contig", "pos", "count_M", "count_F"}; {print $1, $2, $3, $6}' > $OUTPUTcov
