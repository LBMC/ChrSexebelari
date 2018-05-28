~/Programs/DBG2OLC-master/compiled/SparseAssembler g 15 k 51 LD 0 GS 100000 NodeCovTh 1 EdgeCovTh 1 p1 ~/Documents/stage_mbelari/data/sample/MR_350_mapped_to_Mbelari-01182_1.fastq p2 ~/Documents/stage_mbelari/data/sample/MR_350_mapped_to_Mbelari-01182_2.fastq 

~/Programs/DBG2OLC-master/compiled/DBG2OLC k 17 MinLen 500 AdaptiveTh 0.005 KmerCovTh 4 MinOverlap 20 RemoveChimera 0 Contigs Contigs.txt f ~/Documents/stage_mbelari/data/sample/MBELA01182subset.fastq

seqtk seq -a ~/Documents/stage_mbelari/data/sample/MBELA01182subset.fastq > MBELA01182subset.fasta

cat Contigs.txt MBELA01182subset.fasta > ctg_pb.fasta

ulimit -n unlimited

sh ./split_and_run_sparc.sh ~/Documents/stage_mbelari/results/hybrid_test/backbone_raw.fasta ~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC_Consensus_info.txt ~/Documents/stage_mbelari/results/hybrid_test/ctg_pb.fasta ~/Documents/stage_mbelari/results/hybrid_test/consensus_dir 2 > ~/Documents/stage_mbelari/results/hybrid_test/cns_log.txt

sh ./split_and_run_sparc.sh ~/Documents/stage_mbelari/results/hybrid_test/backbone_raw.fasta ~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC_Consensus_info.txt ~/Documents/stage_mbelari/results/hybrid_test/ctg_pb.fasta ~/Documents/stage_mbelari/results/hybrid_test/consensus_dir 2 > ~/Documents/stage_mbelari/results/hybrid_test/cns_log.txt

