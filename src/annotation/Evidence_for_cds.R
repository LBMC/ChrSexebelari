library(seqinr)
source('~/Documents/stage_mbelari/src/Functions.R')

cds <- read.fasta('~/Documents/stage_mbelari/results/annotation/BRAKER/soft_masked/augustus.hints.codingseq')
cds_names <- getName(cds)
cds_seqs <- getSequence(cds)
tmp_names <- sub('Backbone_[0-9]+\\.','',cds_names)

score <- readLines('~/Documents/stage_mbelari/results/annotation/BRAKER/soft_masked/evidence_cds.txt')
score <- sub('# % of transcript supported by hints \\(any source\\): ','',score)
score <- as.numeric(score)
names(score) <- tmp_names
tmp_names <- sub('\\.t[0-9]+','',tmp_names)

matrix_scores <- cbind(tmp_names, score)

genes.absent.females <- readLines('~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-24-Genes_absent_in_females.txt')
genes.poor.females <- readLines('~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-17-Genes_poor_in_females.txt')

#genes.fem.abs.scores <- score[match(genes.poor.females,matrix_scores[,1])]
#new_list <- matrix_scores[names(genes.fem.abs.scores[genes.fem.abs.scores!=0]),]

#genes.fem.absent.scores <- score[match(genes.absent.females,matrix_scores[,1])]
#new_list_absent <- matrix_scores[names(genes.fem.absent.scores[genes.fem.abs.scores!=0]),]

#write.table(new_list[,1],'2018-07-24-Genes_poor_in_females_with_transcript_evidence.txt', quote=F, col.names=F,row.names=F, sep='\t')

getCdsH(new_list[,1],'2018-07-24-Genes_poor_in_females_with_transcript_evidence.fasta')
getCdsH(genes.absent.females,'2018-07-24-Genes_absent_in_females.fasta')

#system("bash ~/Documents/stage_mbelari/src/13genes/Blast-vs-draft-hybrid-genome.sh")

blast_results <- read.table('~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-18-BLAST_Genes_poor_in_females_with_transcript_evidence.txt', head=F, sep='\t')
blast_results <- read.table('~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-24-BLAST_Genes_absent_in_females.txt', head=F, sep='\t')

cov.contig <- read.table('~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-FC_normalized_coverage_at_contig.txt', sep='\t', head=T,row.names=1)

count.hits.contig <- table(blast_results[,1],blast_results[,2])
count.hits.contig.1based <- table(blast_results[,1],blast_results[,2])
count.hits.contig.1based[count.hits.contig.1based>1]<-1
Backbones <- sort(apply(count.hits.contig.1based,2,sum), decreasing=T)


#Create a bed file from the blast hits positions
bed.to.intersect <- blast_results[,c(2,9,10)]
sense <- as.numeric(bed.to.intersect[,2] < bed.to.intersect[,3])
sense[sense==1]<- '+'
sense[sense==0]<- '-'
bed.to.intersect[sense=='-',2] <- blast_results[sense=='-',10]
bed.to.intersect[sense=='-',3] <- blast_results[sense=='-',9]
hits <- sub('Backbone_[0-9]+\\.','',blast_results[,1])
scores <- blast_results$V3 / 100
bed.to.intersect <- data.frame(bed.to.intersect,hits,scores,sense)
bed.to.intersect <- bed.to.intersect[order(bed.to.intersect[,2]),]
bed.to.intersect <- bed.to.intersect[order(bed.to.intersect[,1]),]

write.table(bed.to.intersect,'2018-07-24-Blast_female_absent_genes_hits.bed', quote=F, col.names=F,row.names=F, sep='\t')

system("bash ~/Documents/stage_mbelari/src/coverage_analysis/get_counts_per_blast_hit.sh")

count.blast.hits.fem <- read.table('~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-24-Blast_female_absent_genes_hits_MRDR5_counts.txt', sep='\t', head=F)
count.blast.hits.mal <- read.table('~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-24-Blast_female_absent_genes_hits_MRDR6_counts.txt', sep='\t', head=F)

count.blast.hits.fem.norm <- count.blast.hits.fem
count.blast.hits.mal.norm <- count.blast.hits.mal

count.sum.fem <- c()
count.sum.mal <- c()

for (n in 1:length(Backbones)){
	contig=names(Backbones)[n]
	tmp_fem <- count.blast.hits.fem[,1]==contig
	tmp_mal <- count.blast.hits.mal[,1]==contig
	count.blast.hits.fem.norm[tmp_fem,7] <- (count.blast.hits.fem[tmp_fem,7])*((cov.contig[contig,'counts.norm.female']+1)/(cov.contig[contig,'counts.raw.female']+1))
	count.blast.hits.mal.norm[tmp_mal,7] <- (count.blast.hits.mal[tmp_mal,7])*((cov.contig[contig,'counts.norm.male']+1)/(cov.contig[contig,'counts.raw.male']+1))
	count.fem <- tapply(count.blast.hits.fem.norm[tmp_fem,7],as.character(count.blast.hits.fem.norm[tmp_fem,4]),sum)
	count.mal <- tapply(count.blast.hits.mal.norm[tmp_mal,7],as.character(count.blast.hits.mal.norm[tmp_mal,4]),sum)
	tmp_count_fem <- cbind(rep(contig,length(count.fem)),names(count.fem),count.fem)
	tmp_count_mal <- cbind(rep(contig,length(count.mal)),names(count.mal),count.mal)
	row.names(tmp_count_fem)<-NULL
	row.names(tmp_count_mal)<-NULL
	count.sum.fem <- rbind(count.sum.fem,tmp_count_fem)
	count.sum.mal <- rbind(count.sum.mal,tmp_count_mal)
	print(paste(nrow(count.sum.fem),nrow(count.sum.mal)))
}


#write.table(count.sum.fem,'2018-07-24-normalized-reads-blast-hits-female.txt', sep='\t',row.names=F, col.names=F, quote=F)
#write.table(count.sum.mal,'2018-07-24-normalized-reads-blast-hits-male.txt', sep='\t',row.names=F, col.names=F, quote=F)


#count.sum.fem <- read.table('~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-24-normalized-reads-blast-hits-female.txt', sep='\t', head=F)
#count.sum.mal <- read.table('~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-24-normalized-reads-blast-hits-male.txt', sep='\t', head=F)

count.sum.log2fc <- data.frame(count.sum.fem,count.sum.mal$V3)


fem.reads <- table(count.sum.fem[,1],count.sum.fem[,2])
mal.reads <- table(count.sum.mal[,1],count.sum.mal[,2])
fem.reads <- as.matrix(fem.reads, 'numeric')
mal.reads <- as.matrix(mal.reads, 'numeric')

for (x in 1:nrow(count.sum.fem)){
fem.reads[count.sum.fem[x,1],count.sum.fem[x,2]] <- as.numeric(count.sum.fem[x,3])
mal.reads[count.sum.mal[x,1],count.sum.mal[x,2]] <- as.numeric(count.sum.mal[x,3])
}

log2fc <- log2((fem.reads+1)/(mal.reads+1))
log2fc[fem.reads==0&mal.reads==0]<-NA
#write.table(log2fc,'2018-07-23-log2fc-normalized-reads-blast-hits-female-vs-male.txt', sep='\t',row.names=F, col.names=F, quote=F)
#log2fc<- read.table('~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-23-log2fc-normalized-reads-blast-hits-female-vs-male.txt', sep='\t', head=F)
log2fc[mal.reads<5&fem.reads<10] <- NA
log2fc <- as.matrix(log2fc)
rownames(log2fc) <- sub('Backbone_','B_',rownames(log2fc))


library(gplots)
library(RColorBrewer)
my_col <- colorRampPalette(c("green","black","red"))(256)
heatmap.2(log2fc, col= my_col, trace='none', Rowv=FALSE, Colv=FALSE,keysize=1,key.title=NA, key.xlab='Log2FC', cexRow=0.7, cexCol=1.2, margins=c(5.7,7) , notecol='white', density.info='none')

#cleaner heatmap

clean <- apply(apply(log2fc,2,function(x)is.na(x)),2,sum)
log2fc_clean <- log2fc[,names(which(clean!=100))]
#write.table(log2fc_clean,'2018-07-24-log2fc_clean-normalized-reads-blast-hits-female-vs-male.txt', sep='\t',row.names=T, col.names=T, quote=F)
log2fc_clean<- read.table('~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-24-log2fc_clean-normalized-reads-blast-hits-female-vs-male.txt', sep='\t', head=T, row.names=1)
log2fc_clean <- as.matrix(log2fc_clean)
heatmap.2(log2fc_clean, col= my_col, trace='none', Rowv=FALSE, Colv=FALSE,keysize=1,key.title=NA, key.xlab='Log2FC', cexRow=0.7, cexCol=1.2, margins=c(5.7,7) , notecol='white', density.info='none')
heatmap.2(t(log2fc_clean), col= my_col, trace='none', Rowv=FALSE, Colv=FALSE,keysize=1,key.title=NA, key.xlab='Log2FC', cexRow=0.8, cexCol=0.75, margins=c(5.7,7) , notecol='white', density.info='none')



