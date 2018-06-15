library(seqinr)

#genome contigs
genome <- read.fasta('~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta')
genome_names <- getName(genome)
genome_seqs <- getSequence(genome)

table <- read.csv('Genes_in_scaffolds.csv')
male_backbones <- names(table)[which(table[13,]=='m')]




