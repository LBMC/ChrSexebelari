#!/bin/bash
### shell du job
#$ -S /bin/bash
### nom du job (a changer)
#$ -N bamcompare
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
module load ~/privatemodules
ml deepTools/3.0.2

### execution du programme
###EXECDIR=${HOME}/Formations/Sequentiel
###${EXECDIR}/SommeVecVecSEQ.exe < Monfichierdedata > monfichierresultat

bam1=/scratch/lestrada/stage_mbelari/results/13genes/2018-06-12-MRDR5_vs_Hybrid_assembly.sorted.filtered.bam
bam2=/scratch/lestrada/stage_mbelari/results/13genes/2018-06-12-MRDR6_vs_Hybrid_assembly.sorted.filtered.bam

bamCompare -b1 ${bam1} -b2 ${bam2} -o log2ratio.bed -of bedgraph --ignoreDuplicates



# fin```

