library(seqinr)
genome <- read.fasta('~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta')
genome_names <- getName(genome)
genome_seqs <- getSequence(genome)

genomesizes <- cbind(genome_names,lengths(genome_seqs))

write.table(genomesizes,'2018-06-21-Mbelari_hybrid_genome_sizes.txt', sep='\t', quote=FALSE, col.names=FALSE,row.names=FALSE)
