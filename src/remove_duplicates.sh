#$ -S /bin/bash
### nom du job:
#$ -N remove_duplicates
### file d'attente:
#$ -q E5-2670deb128*
### parallel environnement & nslots
#$ -pe openmp16 16
### charger l'environnement utilisateur pour SGE
#$ -cwd
### exporte les variables d'environnement sur les noeuds d'exécution
#$ -V
### change logs folder
#$ -o /home/cburny/logs
#$ -e /home/cburny/logs

# initialiser environnement Module
source /usr/share/modules/init/bash
module use /applis/PSMN/Modules
module load Base/psmn
module load SAMtools/1.3.1
module list
umask 002

INPUT_BAM="/scratch/cburny/Output_Mbelari/2017_09_12_MRDR6_trim_Mbelari_mapped.bam"
OUTPUT="/scratch/cburny/Output_Mbelari/2017_09_13_MRDR6_trim_Mbelari_mapped_sort"
OUTPUT_BAM="/scratch/cburny/Output_Mbelari/2017_09_13_MRDR6_trim_Mbelari_mapped_sort.bam"
OUTPUTrmdup="/scratch/cburny/Output_Mbelari/2017_09_13_MRDR6_trim_Mbelari_mapped_sort_rmdup.bam"
REPORTrmdup="/scratch/cburny/Output_Mbelari/2017_09_13_MRDR6_trim_Mbelari_mapped_sort_rmdup_flagstat.txt"

##### Sort and index
samtools sort $INPUT_BAM $OUTPUT &&\
samtools index $OUTPUT_BAM &&\

##### remove of duplicates
samtools rmdup $OUTPUT_BAM $OUTPUTrmdup &&\

samtools flagstat $OUTPUTrmdup > $REPORTrmdup
