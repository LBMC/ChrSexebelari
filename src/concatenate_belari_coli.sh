# Concatenate E. coli and M. belari fasta. This was done on uncompressed fasta files
# done on local computer 

cat data/ReferenceGenomes/GCA_000176815.1_ASM17681v1_genomic.fna data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2.scaffolds.fa.repeatmasker.masked | gzip > data/ReferenceGenomes/combined_belari_coli.fasta.gz
