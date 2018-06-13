makeblastdb -in ~/Documents/stage_mbelari/results/hybrid_test/2018-06-01-DGB2OLC/final_assembly.fasta -parse_seqids -dbtype nucl

blastn -query ~/Documents/stage_mbelari/results/13genes/version2/2018-06-11-13genes_genomic_region.fasta \
-db ~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-07/final_assembly_min500_001.fasta -out 2018-06-11-13genes_in_hybrid_genome_v3.txt -outfmt 6

blastn -query ~/Documents/stage_mbelari/results/13genes/version2/2018-06-11-13genes_genomic_region.fasta -db ~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-07/final_assembly_min500_001.fasta -out 2018-06-11-13genes_in_hybrid_genome_v3-align.txt

blastn -query ~/Documents/stage_mbelari/results/13genes/version2/2018-06-11-Contigs_test.fasta \
-db ~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta \
-out 2018-06-11-males_contigs_in_hybrid_genome-align.txt

blastn -query ~/Documents/stage_mbelari/results/13genes/version2/2018-06-11-Contigs_test.fasta \
-db ~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta \
-out 2018-06-11-males_contigs_in_hybrid_genome.txt -outfmt 6