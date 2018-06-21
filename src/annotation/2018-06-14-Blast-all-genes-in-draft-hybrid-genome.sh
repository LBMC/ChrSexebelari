#makeblastdb -in ~/Documents/stage_mbelari/results/hybrid_test/2018-06-07-DGB2OLC/final_assembly_min500_001.fasta -parse_seqids -dbtype nucl

blastn -query ~/Documents/stage_mbelari/data/raw/Mesorhabditis_belari_JU2817_v2.cds.fa -db ~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-07/final_assembly_min500_001.fasta -out 2018-06-14-blast_all_genes_vs_hybrid.txt -outfmt 6

blastn -query ~/Documents/stage_mbelari/data/raw/Mesorhabditis_belari_JU2817_v2.cds.fa -db ~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-07/final_assembly_min500_001.fasta -out 2018-06-14-blast_all_genes_vs_hybrid_align.txt
