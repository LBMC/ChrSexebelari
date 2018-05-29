#!/usr/bin/env bash

illumina_fq_1=$1
illumina_fq_2=$2
nanopore_reads_fq=$3

chunk=`basename ${nanopore_reads_fq} .fastq`

seqtk seq -a ${nanopore_reads_fq} > ${chunk}.fasta

nanopore_reads_fa= ${chunk}.fasta

~/Programs/DBG2OLC-master/compiled/SparseAssembler g 15 k 51 LD 0 GS 100000 NodeCovTh 2 EdgeCovTh 1 p1 ${illumina_fq_1} p2 ${illumina_fq_2}

#Nanopore files and Contigs.txt file from previous step

~/Programs/DBG2OLC-master/compiled/DBG2OLC k 17 MinLen 250 AdaptiveTh 0.001 KmerCovTh 3 MinOverlap 15 RemoveChimera 1 Contigs Contigs.txt f ${nanopore_reads_fq}

#Convert fastq nanopore reads to fasta

cat Contigs.txt ${nanopore_reads_fa} > ctg_pb.fasta

#ulimit -n unlimited

#Consensus step

sh ~/Programs/DBG2OLC-master/utility/split_and_run_sparc_test.sh backbone_raw.fasta DBG2OLC_Consensus_info.txt ctg_pb.fasta consensus_dir 5 > cns_log.txt
