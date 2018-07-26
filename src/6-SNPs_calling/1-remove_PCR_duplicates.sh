INPUTbam=$1
OUTPUTbam=$2

samtools rmdup ${INPUTbam} ${OUTPUTbam}
