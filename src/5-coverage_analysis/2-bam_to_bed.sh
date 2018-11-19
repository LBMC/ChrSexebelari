#The bam files were filtered by quality

INPUTcov=results/coverage_analysis/2018-06-12-MRDR6_vs_Hybrid_assembly.sorted.filtered.bam
GENOMEsizes=results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt
OUTPUTcov=results/coverage_analysis/2018_07_25_MRDR6_vs_Hybrid_assembly.sorted.filtered.bed

bedtools genomecov -d -ibam $INPUTcov -g $GENOMEsizes > $OUTPUTcov

