library(yarrr)
library(seqinr)
cds <- read.fasta('~/Documents/stage_mbelari/data/raw/Mesorhabditis_belari_JU2817_v2.cds.fa')
cds_names <- getName(cds)
cds_seqs <- getSequence(cds)

cdsh <- read.fasta('~/Documents/stage_mbelari/results/annotation/BRAKER/soft_masked/augustus.hints.codingseq')
cdsh_names <- getName(cdsh)
cdsh_seqs <- getSequence(cdsh)
sumh <- summary(lengths(cdsh_seqs))[c(1,3,4,6)]
sumi <- summary(lengths(cds_seqs))[c(1,3,4,6)]

colors <- c(transparent(orig.col = "firebrick1", trans.val = 0.5, maxColorValue = 255),transparent(orig.col = "forestgreen", trans.val = 0.5, maxColorValue = 255))
hist(lengths(cdsh_seqs), breaks=100, col=colors[1], border='white',
main=NULL, xlab='cds length (bp)')
hist(lengths(cds_seqs), breaks=100, col=colors[2], border='white',add=T)
legend('topright',c('Annotation in hybrid assembly genome', 'Annotation in Illumina assembly'), fill=colors, border='white', cex=0.8)
legend('right',c(paste0('Max length: ',summary(lengths(cdsh_seqs))[6],'bp'),paste0('Max length: ',summary(lengths(cds_seqs))[6],'bp')), bty='n', cex=0.8, fill=colors, border='white')


boxplot(lengths(cdsh_seqs), lengths(cds_seqs))




