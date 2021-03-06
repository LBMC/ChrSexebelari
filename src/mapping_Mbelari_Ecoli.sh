#$ -S /bin/bash
### nom du job:
#$ -N mapping_belari_bowtie2
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
module load SAMtools/1.3.1
module list

INDEX="/scratch/cburny/Ref_Mbelari/Ref_Mbelari_Ecoli"
OPTION="-q --very-sensitive -I 100 -X 900"
READS1="/scratch/cburny/Input_ChrSexe/2017_09_11_MRDR5_trim_R1.fastq.gz"
READS2="/scratch/cburny/Input_ChrSexe/2017_09_11_MRDR5_trim_R2.fastq.gz"
OUTPUT="/scratch/cburny/Output_Mbelari/MRDR5_trim_Mbelari_Ecoli.bam"
REPORT_MAPPING="/scratch/cburny/Output_Mbelari/MRDR5_trim_Mbelari_Ecoli_flagstat.txt"
REPORT_BOWTIE2="/scratch/cburny/Output_Mbelari/MRDR5_trim_Mbelari_Ecoli_bowtie2.txt"

# Mapping and BAM output
bowtie2 $OPTION -x $INDEX -1 $READS1 -2 $READS2 2> $REPORT_BOWTIE2 | samtools view -Sb - > $OUTPUT &&\

# Mapping report
samtools flagstat $OUTPUT > $REPORT_MAPPING
