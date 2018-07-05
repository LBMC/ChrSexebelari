library(preprocessCore)
library(yarrr)
fem <- read.table('~/Documents/stage_mbelari/results/coverage_analysis/2018-06-28-g17460_contigs_subbed_5_2.bed', sep='\t', head=F)
male <- read.table('~/Documents/stage_mbelari/results/coverage_analysis/2018-06-28-g17460_contigs_subbed_6_2.bed', sep='\t', head=F)

counts.contigs <- read.table('~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-FC_normalized_coverage_at_contig.txt', sep='\t',head=T,row.names=1)

contigs2graph <- c('1437','207','495','6','61','713','912','11')
contigs2graph <- paste0('Backbone_',contigs2graph)


for (n in 1:length(contigs2graph)){
	contig <- contigs2graph[n]
	tempf <- fem[fem[,1]==contig,]
	tempm <- male[male[,1]==contig,]
	norm.male <- tempm[,3]*(counts.contigs[contig,'counts.norm.male']/counts.contigs[contig,'counts.raw.male'])
	norm.female <- tempf[,3]*(counts.contigs[contig,'counts.norm.female']/counts.contigs[contig,'counts.raw.female'])
	#png(paste0(contig,'_plot.png'))
	plot(norm.male[1:8000], type='l', col=transparent(orig.col = "darkgreen", trans.val = 0.3, maxColorValue = 255), ylim=c(0, 200), xlab='position', ylab='normalized reads', main=contig, lwd=2)
#61	lines(norm.female[503254:511254], type='l', col=transparent(orig.col = 'black', trans.val = 0.3, maxColorValue = 255), lwd=2)
#61	arrows(x0=3771,x1=3771+457,y0=50,y1=50, code=3, angle=90, length=0.05)
#61	text(3771+(457/2),60, labels='g17460', cex=1.5)

#912	lines(norm.female[1:8000], type='l', col=transparent(orig.col = 'black', trans.val = 0.3, maxColorValue = 255), lwd=2)
#912	arrows(x0=3516,x1=3971,y0=150,y1=150, code=3, angle=90, length=0.05)
#912	text(3516+(457/2),155, labels='g17460', cex=1.5)


#495	lines(norm.female[5001:13001], type='l', col=transparent(orig.col = 'black', trans.val = 0.3, maxColorValue = 255), lwd=2)
#495	arrows(x0=3771,x1=3771+457,y0=165,y1=165, code=3, angle=90, length=0.05)
#495	text(3771+(457/2),170, labels='g17460', cex=1.5)

#RNAP	lines(norm.female[8596:16596], type='l', col=transparent(orig.col = 'black', trans.val = 0.3, maxColorValue = 255), lwd=2)
#RNAP	arrows(x0=3501,x1=3501+998,y0=50,y1=50, code=3, angle=90, length=0.05)

#RNAP	text(3501+(998/2),55, labels='RNAP2', cex=1.5)


	#dev.off()

}

