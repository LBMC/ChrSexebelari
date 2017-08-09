#$ -S /bin/bash
### nom du job:
#$ -N mapping_filter_belari_bowtie2
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


INDEX="/scratch/cburny/Input_Mbelari/Mbelari_SOFT"
INPUT_BAM="/scratch/cburny/Output_Mbelari/2017_08_08_MRDR5_trim_Mbelari_SOFT.bam"
OPTION1="-b -F 300"
OPTION2="-b -q 30"
OUTPUT1="/scratch/cburny/Output_Mbelari/2017_08_08_MRDR5_trim_Mbelari_SOFT_mapped_uni.bam"
OUTPUT2="/scratch/cburny/Output_Mbelari/2017_08_08_MRDR5_trim_Mbelari_SOFT_mapped_uni_qual.bam"

samtools view $OPTION1 $INPUT_BAM -o $OUTPUT1
samtools view $OPTION2 $OUTPUT1 -o $OUTPUT2


##### Run samtools on PSMN
OUTPUT1b="/scratch/cburny/Output_Mbelari/2017_08_08_MRDR5_trim_Mbelari_SOFT_mapped_uni_sort.bam"
OUTPUT2b="/scratch/cburny/Output_Mbelari/2017_08_08_MRDR5_trim_Mbelari_SOFT_mapped_uni_qual_sort.bam"

samtools sort -o $OUTPUT1b $OUTPUT1
samtools index $OUTPUT1b

samtools sort -o $OUTPUT2b $OUTPUT2
samtools index $OUTPUT2b

