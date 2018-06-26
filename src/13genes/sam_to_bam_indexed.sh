INPUTsam=~/Documents/stage_mbelari/results/13genes/version2/2016-06-22-Blast-cds-top-50.sam

BAM=`basename ${INPUTsam} .sam`

samtools view -Sb ${INPUTsam} > ${BAM}.bam

#samtools sort ${BAM}.bam -o ${BAM}.sorted.bam

#samtools index ${BAM}.sorted.bam
