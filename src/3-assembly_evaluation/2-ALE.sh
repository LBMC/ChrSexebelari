INPUTbam=~/Documents/stage_mbelari/results/assembly_comparisons/2018-07-02-Illumina_vs_hybrid_improved_chim1_subset.sorted.bam

ASSEMBLY=final_assembly.fasta


OUTPUT=hybrid_assembly

~/Programs/ALE-master/src/ALE ${INPUTbam} ${ASSEMBLY} ALE_analysis_${OUTPUT}.ale


