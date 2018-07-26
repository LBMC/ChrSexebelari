library(seqinr)
contig <- read.fasta('Redundans_longest_contig_8Mb.fasta')
contig_seq <- getSequence(contig)
telomericseq <- gregexpr('ttaggc',contigc)

