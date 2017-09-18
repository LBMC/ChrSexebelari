#!/bin/bash
### shell du job:
#$ -S /bin/bash
### nom du job:
#$ -N PICARD_dict
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
REF="/scratch/cburny/Ref_Mbelari/Mesorhabditis_belari_JU2817_v2_scaffolds_repeatmasker_masked.fa"
REFOUT="/scratch/cburny/Ref_Mbelari/Mesorhabditis_belari_JU2817_v2_scaffolds_repeatmasker_masked.dict"

java -jar $PICARDHOME/picard.jar CreateSequenceDictionary R=$REF O= $REFOUT

echo "termine"
