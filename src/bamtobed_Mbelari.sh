### nom du job:
#$ -N bam_to_bed_belari
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

INPUTcov="/scratch/cburny/Output_Mbelari/2017_09_13_MRDR5_trim_Mbelari_mapped_sort.bam"
GENOMEsizes="/scratch/cburny/Ref_Mbelari/2017_09_13_Mbelari.sizes.genome"
OUTPUTcov="/scratch/cburny/Output_Mbelari/2017_09_13_MRDR5_trim_Mbelari_mapped_sort.bed"

bedtools genomecov -d -ibam $INPUTcov -g $GENOMEsizes > $OUTPUTcov
