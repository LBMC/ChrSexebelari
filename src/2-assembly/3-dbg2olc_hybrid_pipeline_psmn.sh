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

illumina_fq_1_1=/scratch/lestrada/stage_mbelari/data/raw/MR_350_clean_1.fastq
illumina_fq_1_2=/scratch/lestrada/stage_mbelari/data/raw/MR_350_clean_2.fastq 
illumina_fq_2_1=/scratch/lestrada/stage_mbelari/data/raw/MR_550_clean_1.fastq
illumina_fq_2_2=/scratch/lestrada/stage_mbelari/data/raw/MR_550_clean_2.fastq
genome_size=200000000
nanopore_reads_fa=/scratch/lestrada/stage_mbelari/data/2018-05-24-Nanopore_all_reads_best_qual_wo_ecoli.fasta

/scratch/lestrada/stage_mbelari/bin/DBG2OLC-master/compiled/SparseAssembler g 15 k 51 LD 0 GS ${genome_size} NodeCovTh 2 EdgeCovTh 1 p1 ${illumina_fq_1_1} p2 ${illumina_fq_1_2} p1 ${illumina_fq_2_1} p2 ${illumina_fq_2_2}


#Nanopore files and Contigs.txt file from previous step

/scratch/lestrada/stage_mbelari/bin/DBG2OLC-master/compiled/DBG2OLC k 17 MinLen 250 AdaptiveTh 0.001 KmerCovTh 3 MinOverlap 15 RemoveChimera 1 Contigs Contigs.txt f ${nanopore_reads_fa}

cat Contigs.txt ${nanopore_reads_fa} > ctg_pb.fasta

#ulimit -n unlimited

#Consensus step

sh /scratch/lestrada/stage_mbelari/bin/DBG2OLC-master/utility/split_and_run_sparc_test.sh backbone_raw.fasta DBG2OLC_Consensus_info.txt ctg_pb.fasta consensus_dir 2 > cns_log.txt


# fin```
