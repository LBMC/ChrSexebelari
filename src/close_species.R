mal <- read.table('2018-07-11-MRDR3_vs_MRDR4contigs.sorted.bam.counts', sep='\t')
fem <- read.table('2018-07-11-MRDR4_vs_MRDR3contigs.sorted.bam.counts', sep='\t')
fem_tmp<-fem[fem[,2]>=1000,]
mal_tmp<-mal[mal[,2]>=1000,]

hist(sort(mal_tmp[,3])[1:(nrow(mal_tmp)-1000)], breaks=200)
abline(v=30, col='darkred')
hist(sort(fem_tmp[,3])[1:(nrow(fem_tmp)-1000)], breaks=200)
abline(v=30, col='darkred')


