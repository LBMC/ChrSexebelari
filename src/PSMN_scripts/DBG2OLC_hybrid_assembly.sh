#ssembly_test_5.correctedReads.fasta!/bin/bash
### shell du job
#$ -S /bin/bash
### nom du job (a changer)
#$ -N DBG2OLC_improvement
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

### execution du programme
###EXECDIR=${HOME}/Formations/Sequentiel
###${EXECDIR}/SommeVecVecSEQ.exe < Monfichierdedata > monfichierresultat

llumina_fq_1_1=
illumina_fq_1_2=
illumina_fq_2_1=
illumina_fq_2_2=
genome_size=162
nanopore_reads_fa=

/scratch/lestrada/stage_mbelari/src/2018-06-26-DBG2OLC-improve-assembly.parameters


# fin```

