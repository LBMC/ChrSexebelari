INPUTbam=~/Documents/stage_mbelari/results/13genes/2018-06-12-MRDR5_vs_Hybrid_assembly.sorted.filtered.bam
OUTPUT=~/Documents/stage_mbelari/results/coverage_analysis/2018_06_21_MRDR5_vs_Hybrid_assembly.sorted.filtered.counts

samtools idxstats ${INPUTbam} > ${OUTPUT}
