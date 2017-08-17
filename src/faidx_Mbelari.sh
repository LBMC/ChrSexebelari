#$ -S /bin/bash
### nom du job:
#$ -N faidx_belari
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


REF_fa="/scratch/cburny/Input_Mbelari/mesorhabditis_belari_assembly.2.0.SOFTMASKED.fa"

samtools faidx $REF_fa
