makeblastdb -in ~/Documents/stage_mbelari/results/hybrid_test/2018-06-01-DGB2OLC/final_assembly.fasta -parse_seqids -dbtype nucl

blastn -query ~/Documents/stage_mbelari/data/raw/2018-05-18-Complete_13_genes.fasta -db ~/Documents/stage_mbelari/results/hybrid_test/2018-06-01-DGB2OLC/final_assembly.fasta -out 2018-06-01-13genes_in_hybrid_genome.txt -outfmt 6

blastn -query ~/Documents/stage_mbelari/data/raw/2018-05-18-Complete_13_genes.fasta -db ~/Documents/stage_mbelari/results/hybrid_test/2018-06-01-DGB2OLC/final_assembly.fasta -out 2018-06-01-13genes_in_hybrid_genome-align.txt
