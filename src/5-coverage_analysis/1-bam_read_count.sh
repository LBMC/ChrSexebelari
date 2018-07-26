INPUTbam=results/coverage_analysis/2018-06-12-MRDR6_vs_Hybrid_assembly.sorted.filtered.bam
OUTPUT=results/coverage_analysis/2018-06-12-MRDR6_vs_Hybrid_assembly.sorted.filtered.counts

samtools idxstats ${INPUTbam} > ${OUTPUT}
