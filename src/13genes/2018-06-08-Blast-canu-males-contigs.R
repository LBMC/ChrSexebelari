library(seqinr)

#blast output
file <- read.table('~/Documents/stage_mbelari/results/sex_chromosome/2018-06-08-BLAST_canu_vs_contigs_male.txt',header=F)
genome <- read.fasta('~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta')
genome_names <- getName(genome)
genome_seqs <- getSequence(genome)
contig_size <- lengths(genome_seqs)
names(contig_size) <- genome_names


file <- data.frame(file,contigs_size[file[,2]])
borders <- file[file[,9]<3|file[,10]<3,]
borders2 <- file[(abs(file[,9]-file[,13])<3)|(abs(file[,10]-file[,13])<3),]
overl <- cbind(table(borders[,1]),table(borders2[,1]))
overl <- overl[overl[,1]>=1&overl[,2]>=1,]
overl_borders <- borders[which(!is.na(match(borders[,1],rownames(overl)))),]
overl_borders2 <- borders2[which(!is.na(match(borders2[,1],rownames(overl)))),]

borders0 <- rbind(overl_borders,overl_borders2)
borders0 <- borders0[order(borders0[,1]),]

file[file[,1]==unique(borders[,1])[1],]

