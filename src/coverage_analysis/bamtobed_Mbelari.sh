#The bam files were filtered by quality

INPUTcov=~/Documents/stage_mbelari/results/13genes/Illumina_reads_vs_hgenome.sorted.filtered.bam
GENOMEsizes=~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt
OUTPUTcov=~/Documents/stage_mbelari/results/coverage_analysis/2018_06_21_Illumina_vs_Hybrid_assembly.sorted.filtered.bed

bedtools genomecov -d -ibam $INPUTcov -g $GENOMEsizes > $OUTPUTcov

