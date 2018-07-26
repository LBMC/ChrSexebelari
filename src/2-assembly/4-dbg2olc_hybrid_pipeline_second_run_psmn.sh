#!/bin/bash
### shell du job
#$ -S /bin/bash
### nom du job (a changer)
#$ -N DBG2OLC_full_test_1
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

#This file are the contigs obtained from the first step of the previous run, the output from 3-dbg2olc_hybrid_pipeline_psmn.sh
Contigs=/scratch/lestrada/stage_mbelari/results/hybrid_test/full_test/test1/Contigs.txt
#This file contains the nanopore clean reads in fasta and the output1 of the dbg2olc pipeline (final_assembly.fasta)
nanopore_reads_fa=/scratch/lestrada/stage_mbelari/results/hybrid_test/full_test/test1/2018-06-27-merged_nanopore_reads_hybrid_assembly.fasta

#AdaptiveTh 0.005 

/scratch/lestrada/stage_mbelari/bin/DBG2OLC-master/compiled/DBG2OLC k 17 MinLen 500 AdaptiveTh 0.005 KmerCovTh 3 MinOverlap 15 RemoveChimera 1 Contigs ${Contigs} f ${nanopore_reads_fa}

cat ${Contigs} ${nanopore_reads_fa} > ctg_pb.fasta

#ulimit -n unlimited

#Consensus step

sh /scratch/lestrada/stage_mbelari/bin/DBG2OLC-master/utility/split_and_run_sparc_test.sh backbone_raw.fasta DBG2OLC_Consensus_info.txt ctg_pb.fasta consensus_dir 2 > cns_log.txt


# fin```
