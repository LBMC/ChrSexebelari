INPUTbam=~/Documents/stage_mbelari/results/13genes/Illumina_reads_vs_hgenome.sorted.bam
OUTPUTfilterbam=~/Documents/stage_mbelari/results/13genes/Illumina_reads_vs_hgenome.sorted.filtered.bam

samtools view -bq 5 ${INPUTbam} > ${OUTPUTfilterbam}
