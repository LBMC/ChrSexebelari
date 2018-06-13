makeblastdb -in ~/Documents/stage_mbelari/results/13genes/2018-06-08-Contigs_males.fasta -parse_seqids -dbtype nucl

blastn -query ~/Documents/stage_mbelari/results/canu/full_assembly/2018-05-24-assembly_test_5.complete_assembly.fasta -db ~/Documents/stage_mbelari/results/13genes/2018-06-08-Contigs_males.fasta -out 2018-06-08-BLAST_canu_vs_contigs_male.txt -outfmt 6

blastn -query ~/Documents/stage_mbelari/results/canu/full_assembly/2018-05-24-assembly_test_5.complete_assembly.fasta -db ~/Documents/stage_mbelari/results/13genes/2018-06-08-Contigs_males.fasta -out 2018-06-08-BLAST_canu_vs_contigs_male-align.txt

