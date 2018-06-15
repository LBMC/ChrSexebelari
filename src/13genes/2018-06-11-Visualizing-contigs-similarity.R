library(seqinr)
file <- read.table('~/Documents/stage_mbelari/results/13genes/version2/2018-06-11-males_contigs_in_hybrid_genome.txt',header=F)

contigs_name <- c(207,912,6,713,1437,11,115,117)
contigs_name <- paste0('Backbone_',contigs_name)

pdf('Contigs_blast.pdf')
for (n in 1:length(contigs_name)){
	x0=file[which(!is.na(match(file[,1],contigs_name[n]))),c(7)]
	x0=x0[-1]
	x1=file[which(!is.na(match(file[,1],contigs_name[n]))),c(8)]
	x1=x1[-1]
	y0=seq(1,length(x0))
	plot(x0,y0,type='n',ylab='number of hits',xlab='contig position', col='blue',main=contigs_name[n])	
	arrows(x0=x0,x1=x1,y0=y0,y1=y0,length=0)
	}
dev.off()