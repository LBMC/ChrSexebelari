#$ -S /bin/bash
### nom du job:
#$ -N extract_unmapped_reads_belari_bowtie2
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
module use /applis/PSMN/Modules
module load Base/psmn
module load SAMtools/1.3.1
module list

INPUT_BAM="/scratch/cburny/Output_Mbelari/2017_09_12_MRDR6_trim_Mbelari_Ecoli.bam"
OPTION1="-b -f 4"
OUTPUT1="/scratch/cburny/Output_Mbelari/2017_10_30_MRDR6_trim_Mbelari_Ecoli_unmapped.bam"
REPORTbelari="/scratch/cburny/Output_Mbelari/2017_10_30_MRDR6_trim_Mbelari_Ecoli_unmapped_flagstat.txt"

##### Filter to keep unmapped reads paired or not
samtools view $OPTION1 $INPUT_BAM -o $OUTPUT1 &&\

samtools flagstat $OUTPUT1 > $REPORTbelari
