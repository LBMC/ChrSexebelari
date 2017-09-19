#$ -S /bin/bash
### nom du job:
#$ -N bam_to_bedpe_belari
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
module load BEDtools/2.24.0
module list

OPTION="-bedpe"
INPUT_BAM="/scratch/cburny/Output_Mbelari/2017_09_13_MRDR5_trim_Mbelari_mapped.bam"
OUTPUT="/scratch/cburny/Output_Mbelari/2017_09_13_MRDR5_trim_Mbelari_mapped_sort_read_names"
INPUTcov="/scratch/cburny/Output_Mbelari/2017_09_13_MRDR5_trim_Mbelari_mapped_sort_read_names.bam"
OUTPUTcov="/scratch/cburny/Output_Mbelari/2017_09_13_MRDR5_trim_Mbelari_mapped_sort_read_names.bedpe"

##### Order by read names 
samtools sort -n $INPUT_BAM $OUTPUT &&\

bedtools bamtobed $OPTION -i $INPUTcov > $OUTPUTcov


