INPUTbam=~/Documents/stage_mbelari/results/sex_chromosome/2018-06-29-redundans_assembly_2_MRDR6.sorted.bam
OUTPUT=~/Documents/stage_mbelari/results/coverage_analysis/2018_06_29_MRDR6_vs_redundans_assembly.sorted.filtered.counts

samtools idxstats ${INPUTbam} > ${OUTPUT}
