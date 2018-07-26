#!/bin/bash
### shell du job
#$ -S /bin/bash
### nom du job (a changer)
#$ -N Redundans_0.8
### file d'attente (a changer)
#$ -q E5-2697Av4deb256
### charger l'environnement utilisateur pour SGE
#$ -cwd
#$ -pe openmp32 4
### exporter les variables d'environnement sur tous les noeuds d'execution
#$ -V
### mails en debut et fin d'execution
#$ -m be
# aller dans le repertoire de travail/soumission
# important, sinon, le programme est lancé depuis ~/
cd ${SGE_O_WORKDIR}
echo ${SGE_O_WORKDIR}
echo ${HOSTNAME}

### configurer l'environnement
source /usr/share/lmod/lmod/init/bash

### execution du programme
###EXECDIR=${HOME}/Formations/Sequentiel
###${EXECDIR}/SommeVecVecSEQ.exe < Monfichierdedata > monfichierresultat
~/Programs/redundans/redundans.py -v -i /scratch/lestrada/stage_mbelari/data/raw/MR_350_clean_1.fastq /scratch/lestrada/stage_mbelari/data/raw/MR_350_clean_2.fastq /scratch/lestrada/stage_mbelari/data/raw/MR_550_clean_1.fastq /scratch/lestrada/stage_mbelari/data/raw/MR_550_clean_2.fastq -l /scratch/lestrada/stage_mbelari/data/2018-05-24-Nanopore_all_reads_best_qual_wo_ecoli.fasta -f /scratch/lestrada/stage_mbelari/results/hybrid_test/full_test/test1/final_assembly.fasta -o test2 --identity 0.8


# fin```
