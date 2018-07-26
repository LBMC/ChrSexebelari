makeblastdb -in ~/Documents/stage_mbelari/data/raw/mbela_assembly_JU2817.fa -parse_seqids -dbtype nucl

blastn -query ~/Documents/stage_mbelari/data/raw/2018-05-18-Complete_13_genes.fasta -db ~/Documents/stage_mbelari/data/raw/mbela_assembly_JU2817.fa -out 2018-06-05-13genes_in_illumina_genome.txt -outfmt 6

blastn -query ~/Documents/stage_mbelari/data/raw/2018-05-18-Complete_13_genes.fasta -db ~/Documents/stage_mbelari/data/raw/mbela_assembly_JU2817.fa -out 2018-06-05-13genes_in_illumina_genome-align.txt
