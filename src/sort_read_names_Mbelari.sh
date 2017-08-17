#$ -S /bin/bash
### nom du job:
#$ -N sort_belari
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

INPUT_BAM="/scratch/cburny/Output_Mbelari/2017_08_08_MRDR5_trim_Mbelari_SOFT_mapped_qual.bam"
OUTPUT="/scratch/cburny/Output_Mbelari/2017_08_08_MRDR5_trim_Mbelari_SOFT_mapped_qual_sort_read_names"

samtools sort -n $INPUT_BAM $OUTPUT
