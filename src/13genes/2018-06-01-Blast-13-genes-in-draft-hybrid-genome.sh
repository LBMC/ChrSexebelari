seqtk seq -a data/raw/2017_09_13_MRDR5_trim_Mbelari_mapped_sort.fastq.gz > results/13genes/2017_09_13_MRDR5_trim_Mbelari_mapped_sort.fasta

seqtk seq -a data/raw/2017_09_13_MRDR6_trim_Mbelari_mapped_sort.fastq.gz > results/13genes/2017_09_13_MRDR6_trim_Mbelari_mapped_sort.fasta

makeblastdb -in results/13genes/2017_09_13_MRDR6_trim_Mbelari_mapped_sort.fasta -parse_seqids -dbtype nucl

makeblastdb -in results/13genes/2017_09_13_MRDR5_trim_Mbelari_mapped_sort.fasta -parse_seqids -dbtype nucl

blastn -query data/raw/2018-05-18-Complete_13_genes.fasta -db 2017_09_13_MRDR6_trim_Mbelari_mapped_sort.fasta -out 2018-05-19-13genes_in_MRDR6.txt -outfmt 6

blastn -query data/raw/2018-05-18-Complete_13_genes.fasta -db 2017_09_13_MRDR5_trim_Mbelari_mapped_sort.fasta -out 2018-05-19-13genes_in_MRDR5.txt -outfmt 6
