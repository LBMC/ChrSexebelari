#!/bin/bash
### shell du job
#$ -S /bin/bash
### nom du job (a changer)
#$ -N illum_vs_reduced
### file d'attente (a changer)
#$ -q monointeldeb128
### charger l'environnement utilisateur pour SGE
#$ -cwd
### exporter les variables d'environnement sur tous les noeuds d'execution
#$ -V
### mails en debut et fin d'execution
#$ -m be

# aller dans le repertoire de travail/soumission
# important, sinon, le programme est lancé depuis ~/
cd ${SGE_O_WORKDIR}

### configurer l'environnement
source /usr/share/lmod/lmod/init/bash
module load GCC/6.4.0/Bowtie2/2.3.3

### execution du programme
###EXECDIR=${HOME}/Formations/Sequentiel
###${EXECDIR}/SommeVecVecSEQ.exe < Monfichierdedata > monfichierresultat

assembly='/scratch/lestrada/stage_mbelari/results/redundans/redundans/test1/scaffolds.reduced.fa'
pe1='/scratch/lestrada/stage_mbelari/data/raw/s_MR_550_clean_1.fastq'
pe2='/scratch/lestrada/stage_mbelari/data/raw/s_MR_550_clean_2.fastq'
index_name='hybrid_reduced'
sam_name='2016-06-26-Illumina_vs_reduced_assembly_default.sam'
bam=`basename ${sam_name} .sam`

bowtie2-build ${assembly} ${index_name}

bowtie2 --very-sensitive-local -x ${index_name} -1 ${pe1} -2 ${pe2} -S ${sam_name}

samtools view -Sb ${sam_name} > ${bam}.bam

samtools sort ${bam}.bam -o ${bam}.sorted.bam

samtools index ${bam}.sorted.bam



# fin```

