#!/bin/bash
### shell du job:
#$ -S /bin/bash
### nom du job:
#$ -N GATK_create_intervals_realign
### file d'attente:
#$ -q E5-2670deb128* 
### parallel environnement & nslots
#$ -pe mpi16_debian 16 
### charger l'environnement utilisateur pour SGE
#$ -cwd
### exporte les variables d'environnement sur les noeuds d'exécution
#$ -V
### change logs folder
#$ -o /home/cburny/logs
#$ -e /home/cburny/logs

# initialiser environnement Module
unset module
source /usr/share/modules/init/bash
module use /applis/PSMN/Modules
module load Base/psmn
module load GATK/3.6 
module list
umask 002

GATKHOME="/applis/PSMN/generic/GATK/3.6"
REF="/scratch/cburny/Ref_Mbelari/Mesorhabditis_belari_JU2817_v2_scaffolds_repeatmasker_masked.fa"
OUTPUT="/scratch/cburny/Output_Mbelari/2017_09_18_MRDR5_Mbelari.realign.intervals"
INDELS="/scratch/cburny/Output_Mbelari/2017_09_15_MRDR5_trim_Mbelari_mapped_sort_rmdup_indels.recode.vcf"
INPUT="/scratch/cburny/Output_Mbelari/2017_09_18_MRDR5_trim_Mbelari_mapped_sort_rmdup.bam"

java -jar $GATKHOME/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF -known $INDELS -I $INPUT -o $OUTPUT

echo "termine"
