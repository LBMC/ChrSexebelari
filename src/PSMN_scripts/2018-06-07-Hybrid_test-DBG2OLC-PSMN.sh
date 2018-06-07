#!/usr/bin/env bash

illumina_fq_1_1=$1
illumina_fq_1_2=$2
illumina_fq_2_1=$3
illumina_fq_2_2=$4
genome_size=$5
nanopore_reads_fa=$6

#/scratch/lestrada/stage_mbelari/bin/DBG2OLC-master/compiled/SparseAssembler g 15 k 63 LD 0 GS ${genome_size} NodeCovTh 2 EdgeCovTh 1 p1 ${illumina_fq_1_1} p2 ${illumina_fq_1_2} p1 ${illumina_fq_2_1} p2 ${illumina_fq_2_2}


#Nanopore files and Contigs.txt file from previous step

/scratch/lestrada/stage_mbelari/bin/DBG2OLC-master/compiled/DBG2OLC k 17 MinLen 500 AdaptiveTh 0.005 KmerCovTh 3 MinOverlap 15 RemoveChimera 1 Contigs Contigs.txt f ${nanopore_reads_fa}

cat Contigs.txt ${nanopore_reads_fa} > ctg_pb.fasta

#ulimit -n unlimited

#Consensus step

sh /scratch/lestrada/stage_mbelari/bin/DBG2OLC-master/utility/split_and_run_sparc_test.sh backbone_raw.fasta DBG2OLC_Consensus_info.txt ctg_pb.fasta consensus_dir 2 > cns_log.txt
