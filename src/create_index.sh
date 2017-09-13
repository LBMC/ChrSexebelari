#!/bin/bash
### nom du job:
#$ -N index_bowtie2
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

REF_FA="/scratch/cburny/Ref_Mbelari/2017_09_08_combined_belari_coli.fasta"
INDEX="/scratch/cburny/Ref_Mbelari/Ref_Mbelari_Ecoli"
OUTPUT_INDEX="/scratch/cburny/Ref_Mbelari/Ref_Mbelari_Ecoli_index_bowtie2.txt"

bowtie2-build $REF_FA $INDEX > $OUTPUT_INDEX
