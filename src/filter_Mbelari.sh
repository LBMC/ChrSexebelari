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
module use /applis/PSMN/Modules
module load Base/psmn
module load Bowtie/2.2.4
module load SAMtools/1.3.1
module list

INDEX="/scratch/cburny/Ref_Mbelari/Ref_Mbelari_Ecoli"
INPUT_BAM="/scratch/cburny/Output_Mbelari/2017_09_12_MRDR6_trim_Mbelari_Ecoli.bam"
OPTION1="-b -F 4"
OUTPUT1="/scratch/cburny/Output_Mbelari/2017_09_12_MRDR6_trim_Mbelari_Ecoli_mapped.bam"
OUTPUT1sort="/scratch/cburny/Output_Mbelari/2017_09_12_MRDR6_trim_Mbelari_Ecoli_mapped_sort"
OUTPUT1sortbam="/scratch/cburny/Output_Mbelari/2017_09_12_MRDR6_trim_Mbelari_Ecoli_mapped_sort.bam"
OUTPUTcoli="/scratch/cburny/Output_Mbelari/2017_09_12_MRDR6_trim_Ecoli_mapped.bam"
OUTPUTbelari="/scratch/cburny/Output_Mbelari/2017_09_12_MRDR6_trim_Mbelari_mapped.bam"
REPORTcoli="/scratch/cburny/Output_Mbelari/2017_09_12_MRDR6_trim_Ecoli_mapped_flagstat.txt"
REPORTbelari="/scratch/cburny/Output_Mbelari/2017_09_12_MRDR6_trim_Mbelari_mapped_flagstat.txt"

##### Filter to keep mapped reads paired or not
samtools view $OPTION1 $INPUT_BAM -o $OUTPUT1 &&\

##### Sort and index
samtools sort $OUTPUT1 $OUTPUT1sort &&\
samtools index $OUTPUT1sortbam &&\

##### Store reads mapped to either E. coli or M. belari
contigs_coli=`awk '{print $0}' /scratch/cburny/Ref_Mbelari/2017_09_12_contigs_short_name_Ecoli.txt` &&\
contigs_belari=`awk '{print $0}' /scratch/cburny/Ref_Mbelari/2017_09_12_contigs_short_name_Mbelari.txt` &&\

samtools view -b $OUTPUT1sortbam $contigs_coli > $OUTPUTcoli &&\
samtools view -b $OUTPUT1sortbam $contigs_belari > $OUTPUTbelari &&\

samtools flagstat $OUTPUTcoli > $REPORTcoli &&\
samtools flagstat $OUTPUTbelari > $REPORTbelari
