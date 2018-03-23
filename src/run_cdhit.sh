#$ -S /bin/bash
### nom du job:
#$ -N CDHIT
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

source /usr/share/modules/init/bash
module use /home/cburny/privatemodules_contrib/modulefiles/
module load CDHIT/4.6.8

cdhit-est-2d -i /scratch/cburny/CDHIT_compare_fasta/2017_11_06_MRDR5_trim_Mbelari_Ecoli_unmapped.fa -i2 /scratch/cburny/CDHIT_compare_fasta/2017_11_06_MRDR6_trim_Mbelari_Ecoli_unmapped.fa -o /scratch/cburny/CDHIT_compare_fasta/2017_11_21_cdhit.fa -n 9 -M 0
