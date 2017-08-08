#!/bin/bash
### shell du job:
#$ -S /bin/bash
### nom du job:
#$ -N mapping_coli_bowtie2
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
module load Bowtie/2.3.0
module load SAMtools/1.3.1
module list

REF="/home/cburny/Ref_coli/Ecoli_K12.fa.gz"
INDEX="/home/cburny/Ref_coli/Ecoli_K12.index"
READS1="/home/cburny/Clean_Reads_ChrSexe/2017_08_07_MRDR5_trim_R1.fastq.gz"
READS2="/home/cburny/Clean_Reads_ChrSexe/2017_08_07_MRDR5_trim_R2.fastq.gz"
OPTION="-q --very-sensitive --no-discordant"
OUTPUT="/scratch/cburny/Output_Ecoli/2017_08_07_MRDR5_trim.bam"
REPORT_MAPPING="/scratch/cburny/Output_Ecoli/2017_08_07_MRDR5_trim_flagstat.txt"
REPORT_BOWTIE2="/scratch/cburny/Output_Ecoli/2017_08_07_MRDR5_trim_bowtie2.txt"

# Build index
bowtie2-build --threads 16 $REF $INDEX

# Mapping and BAM output
bowtie2 $OPTION -p 16 -x $INDEX -1 $READS1 -2 $READS2 2> $REPORT_BOWTIE2 | samtools view -Sb - > $OUTPUT 

# Mapping report
samtools flagstat $OUTPUT > $REPORT_MAPPING






