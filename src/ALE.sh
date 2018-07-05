INPUTbam=~/Documents/stage_mbelari/results/assembly_comparisons/2018-07-02-Illumina_vs_hybrid_improved_chim1_subset.sorted.bam

ASSEMBLY=~/Documents/stage_mbelari/results/post_assembly/improve_assembly/final_assembly_improved.fasta
#ASSEMBLY=~/Documents/stage_mbelari/results/redundans/redundans_stats_changed/scaffolds.reduced.fa


OUTPUT=hybrid_improved

~/Programs/ALE-master/src/ALE ${INPUTbam} ${ASSEMBLY} ALE_analysis_${OUTPUT}.ale


