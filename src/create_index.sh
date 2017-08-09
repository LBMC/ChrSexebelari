#!/bin/bash
### nom du job:
#$ -N index_coli_bowtie2
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
module load Bowtie/2.2.4
module list

REF_GZ="/scratch/cburny/Input_Ecoli/Ecoli_K12.fa.gz"
REF_FA="/scratch/cburny/Input_Ecoli/Ecoli_K12.fa"
INDEX="/scratch/cburny/Input_Ecoli/Ecoli_K12"
OUTPUT_INDEX="/scratch/cburny/Input_Ecoli/index_bowtie2.txt"

# Build index
gzip -dc $REF_GZ > $REF_FA

bowtie2-build $REF_FA $INDEX > $OUTPUT_INDEX
