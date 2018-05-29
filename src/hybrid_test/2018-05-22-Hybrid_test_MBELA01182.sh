#Files

#Illumina files 1 and 2 for contigs formation


~/Programs/DBG2OLC-master/compiled/SparseAssembler g 15 k 51 LD 0 GS 100000 NodeCovTh 2 EdgeCovTh 1 p1 ~/Documents/stage_mbelari/data/sample/MR_350_mapped_to_Mbelari-01182_1.fastq p2 ~/Documents/stage_mbelari/data/sample/MR_350_mapped_to_Mbelari-01182_2.fastq 

#Nanopore files and Contigs.txt file from previous step

~/Programs/DBG2OLC-master/compiled/DBG2OLC k 17 MinLen 250 AdaptiveTh 0.001 KmerCovTh 3 MinOverlap 15 RemoveChimera 1 Contigs ~/Documents/stage_mbelari/results/hybrid_test/Contigs.txt f ~/Documents/stage_mbelari/data/sample/MBELA01182subset.fastq

#Convert fastq nanopore reads to fasta

#seqtk seq -a ~/Documents/stage_mbelari/data/sample/MBELA01182subset.fastq > ~/Documents/stage_mbelari/results/hybrid_test/MBELA01182subset.fasta


cat ~/Documents/stage_mbelari/results/hybrid_test/Contigs.txt ~/Documents/stage_mbelari/results/hybrid_test/MBELA01182subset.fasta > ~/Documents/stage_mbelari/results/hybrid_test/ctg_pb.fasta

#ulimit -n unlimited

#Consensus step

sh ~/Programs/DBG2OLC-master/utility/split_and_run_sparc_test.sh ~/Documents/stage_mbelari/results/hybrid_test/backbone_raw.fasta ~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC_Consensus_info.txt ~/Documents/stage_mbelari/results/hybrid_test/ctg_pb.fasta ~/Documents/stage_mbelari/results/hybrid_test/consensus_dir 2 > ~/Documents/stage_mbelari/results/hybrid_test/cns_log.txt

#Analysis of assembly quality

~/Programs/quast-4.6.3/quast.py ~/Documents/stage_mbelari/results/hybrid_test/consensus_dir/final_assembly.fasta -R ~/Documents/stage_mbelari/data/sample/MBELA01182_RNAP2.fa -o ~/Documents/stage_mbelari/results/hybrid_test/assembly_analysis/
