#$ -S /bin/bash
### nom du job:
#$ -N call_var_belari
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
module load SAMtools/1.1
module load BCFtools/1.1
module list

##### Run samtools and bcftools

OPTION="-S -ABQ0"

REF="-u -f /scratch/cburny/Input_Mbelari/mesorhabditis_belari_assembly.2.0.SOFTMASKED.fa"

INPUT="/scratch/cburny/Output_Mbelari/2017_08_08_MRDR5_trim_Mbelari_SOFT_mapped_qual_sort_read_names.bam"

OUTPUT="/scratch/cburny/Output_Mbelari/2017_08_08_MRDR5_trim_Mbelari_SOFT_mapped_qual_sort_raw.vcf"

samtools mpileup $OPTION $REF $INPUT | bcftools call -mv > $OUTPUT
