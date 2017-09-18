#!/bin/bash
### shell du job:
#$ -S /bin/bash
### nom du job:
#$ -N PICARD_RGtag
### file d'attente:
#$ -q E5-2670deb128*
### parallel environnement & nslots
#$ -pe mpi16_debian 16 
### charger l'environnement utilisateur pour SGE
#$ -cwd
### exporte les variables d'environnement sur les noeuds d'exécution
#$ -V
### mails en debut et fin d execution
#$ -m be
### change logs folder
#$ -o /home/cburny/logs
#$ -e /home/cburny/logs

# initialiser environnement Module
source /usr/share/modules/init/bash
module use /applis/PSMN/Modules
module load Base/psmn 
module load picard/2.3.0 
module list
umask 002

PICARDHOME="/applis/PSMN/generic/picard/2.3.0"
INPUT="/scratch/cburny/Output_Mbelari/2017_09_18_MRDR6_trim_Mbelari_mapped_sort_rmdup.bam"
OUTPUT="/scratch/cburny/Output_Mbelari/2017_09_18_MRDR6_trim_Mbelari_mapped_rmdup_rg.bam"

java -jar $PICARDHOME/picard.jar AddOrReplaceReadGroups I=$INPUT O=$OUTPUT RGID=HKCVHBBXX.1 RGLB=male RGPL=illumina RGPU=HKCVHBBXX.1 RGSM=MRDR6 VALIDATION_STRINGENCY=LENIENT

echo "termine"
