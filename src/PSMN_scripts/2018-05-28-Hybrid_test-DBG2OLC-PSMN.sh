#!/usr/bin/env bash

illumina_fq_1=$1
illumina_fq_2=$2
genome_size=$3
nanopore_reads_fa=$4

~/stage_mbelari/bin/DBG2OLC-master/compiled/SparseAssembler g 15 k 51 LD 0 GS ${genome_size} NodeCovTh 2 EdgeCovTh 1 p1 ${illumina_fq_1} p2 ${illumina_fq_2}

#mv Contigs.txt Contigs_350.txt

#~/stage_mbelari/bin/DBG2OLC-master/compiled/SparseAssembler g 15 k 51 LD 0 GS ${genome_size} NodeCovTh 2 EdgeCovTh 1 p1 ${illumina_fq_1} p2 ${illumina_fq_2}

#mv Contigs.txt Contigs_550.txt

#sed -i 's/Contig/Contig_350/g' Contigs_350.txt
#sed -i 's/Contig/Contig_550/g' Contigs_550.txt
#cat Contigs_350.txt Contigs_550.txt > Contigs.txt

#Nanopore files and Contigs.txt file from previous step

~/stage_mbelari/bin/DBG2OLC-master/compiled/DBG2OLC k 17 MinLen 250 AdaptiveTh 0.001 KmerCovTh 3 MinOverlap 15 RemoveChimera 1 Contigs Contigs.txt f ${nanopore_reads_fa}

cat Contigs.txt ${nanopore_reads_fa} > ctg_pb.fasta

#ulimit -n unlimited

#Consensus step

sh ~/stage_mbelari/bin/DBG2OLC-master/utility/split_and_run_sparc_test.sh backbone_raw.fasta DBG2OLC_Consensus_info.txt ctg_pb.fasta consensus_dir 5 > cns_log.txt
