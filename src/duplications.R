source('~/Documents/stage_mbelari/src/functions/Functions.R')
library(seqinr)
library(gplots)
library(ggplot2)
library(RColorBrewer)

#get a small version of gff3 augustus file
system("bash grep 'gene' ~/Documents/stage_mbelari/results/annotation/BRAKER/soft_masked/augustus.hints.gff3 > ~/Documents/stage_mbelari/results/duplications/gene_locations.txt")

gene.gff <- read.table('~/Documents/stage_mbelari/results/duplications/gene_locations.txt', head=F, sep='\t')
gene.gff$V9<-sub('ID=','',gene.gff$V9)
gene.gff$V9<-sub('\\;','',gene.gff$V9)

genome <- read.fasta('~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta')
genome_names <- getName(genome)
genome_seqs <- getSequence(genome)

genomic_seq <- list()

for (n in 1:nrow(gene.gff)){
contig <- as.character(gene.gff$V1)[n]
genomic_seq[[n]] <- genome_seqs[genome_names==contig][[1]][gene.gff$V4[n]:gene.gff$V5[n]]
}


#Produce fasta file with whole genomic regions
#write.fasta(genomic_seq,gene.gff$V9,'~/Documents/stage_mbelari/results/duplications/mblelari_complete_genomic_genes.fasta')

counts.genes.bp <- read.csv("2018-07-25-tests_FC_normalized_coverage_at_bp_within_genes.txt", 
  sep = "\t", h = T, stringsAsFactors = T)
test <- apply(counts.genes.bp,2,function(x){is.na(x)})


genes.absent.females <- readLines('~/Documents/stage_mbelari/results/annotation/coverage_analysis/2018-07-24-Genes_absent_in_females.txt')

#add RNAP2
genes.absent.females <- c(genes.absent.females,'g25986')

genomique.regions <- read.fasta('~/Documents/stage_mbelari/results/duplications/mblelari_complete_genomic_genes.fasta')
genreg_seqs <- getSequence(genomique.regions)
genreg_names <- getName(genomique.regions)

#Produce fasta file of genes absent in females with genomic regions
write.fasta(genreg_seqs[match(genes.absent.females,genreg_names)],genes.absent.females,'~/Documents/stage_mbelari/results/duplications/2018-07-26-complete_seq_genes_absent_females.fasta')

#Blast the sequences
system("bash ~/Documents/stage_mbelari/src/functions/blastn.sh")

blast_results <- read.table('2018-07-26-BLAST_Genes_absent_in_females_genomic_seq.txt', head=F, sep='\t')

genes.absent.fem <- read.fasta('~/Documents/stage_mbelari/results/duplications/2018-07-26-complete_seq_genes_absent_females.fasta')
genes.absent.fem.names <- getName(genes.absent.fem)
genes.absent.fem.seqs <- getSequence(genes.absent.fem)

gene.sizes <- lengths(genes.absent.fem)

cover <- c()

for (x in 1:length(gene.sizes)){
gene <- names(gene.sizes)[x]
cover[blast_results$V1==gene] <- (blast_results$V8[blast_results$V1==gene]-blast_results$V7[blast_results$V1==gene]+1)/gene.sizes[x]
}

blast_results <- data.frame(blast_results,cover)
blast_results_filt <- blast_results[cover>0.8,]


cov.contig <- read.table('~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-FC_normalized_coverage_at_contig.txt', sep='\t', head=T,row.names=1)

count.hits.contig <- table(blast_results_filt[,1],blast_results_filt[,2])
count.hits.contig.1based <- table(blast_results_filt[,1],blast_results_filt[,2])
count.hits.contig.1based[count.hits.contig.1based>1]<-1
Backbones <- sort(apply(count.hits.contig.1based,2,sum), decreasing=T)


#Create a bed file from the blast hits positions
bed.to.intersect <- blast_results_filt[,c(2,9,10)]
sense <- as.numeric(bed.to.intersect[,2] < bed.to.intersect[,3])
sense[sense==1]<- '+'
sense[sense==0]<- '-'
bed.to.intersect[sense=='-',2] <- blast_results_filt[sense=='-',10]
bed.to.intersect[sense=='-',3] <- blast_results_filt[sense=='-',9]
hits <- sub('Backbone_[0-9]+\\.','',blast_results_filt[,1])
scores <- blast_results_filt$V3 / 100
cover <- blast_results_filt$cover
bed.to.intersect <- data.frame(bed.to.intersect,hits,scores,sense,cover)
bed.to.intersect <- bed.to.intersect[order(bed.to.intersect[,2]),]
bed.to.intersect <- bed.to.intersect[order(bed.to.intersect[,1]),]

#write.table(bed.to.intersect,'2018-07-26-Blast_female_absent_genes_hits.bed', quote=F, col.names=F,row.names=F, sep='\t')

#intersect the bed file from blast to know the counts in aligned regions
system("bash ~/Documents/stage_mbelari/src/coverage_analysis/get_counts_per_blast_hit.sh")


#Analyze differences between female and male
count.blast.hits.fem <- read.table('~/Documents/stage_mbelari/results/duplications/2018-07-26-Blast_female_absent_genes_hits_MRDR5_counts.txt', sep='\t', head=F)
count.blast.hits.mal <- read.table('~/Documents/stage_mbelari/results/duplications/2018-07-26-Blast_female_absent_genes_hits_MRDR6_counts.txt', sep='\t', head=F)


#Normalize raw counts
count.blast.hits.fem.norm <- count.blast.hits.fem
count.blast.hits.mal.norm <- count.blast.hits.mal

count.sum.fem <- c()
count.sum.mal <- c()

for (n in 1:length(Backbones)){
	contig=names(Backbones)[n]
	tmp_fem <- count.blast.hits.fem[,1]==contig
	tmp_mal <- count.blast.hits.mal[,1]==contig
	count.blast.hits.fem.norm[tmp_fem,8] <- (count.blast.hits.fem[tmp_fem,8])*((cov.contig[contig,'counts.norm.female']+1)/(cov.contig[contig,'counts.raw.female']+1))
	count.blast.hits.mal.norm[tmp_mal,8] <- (count.blast.hits.mal[tmp_mal,8])*((cov.contig[contig,'counts.norm.male']+1)/(cov.contig[contig,'counts.raw.male']+1))
	#count.fem <- tapply(count.blast.hits.fem.norm[tmp_fem,8],as.character(count.blast.hits.fem.norm[tmp_fem,4]),sum)
	#count.mal <- tapply(count.blast.hits.mal.norm[tmp_mal,8],as.character(count.blast.hits.mal.norm[tmp_mal,4]),sum)
	#tmp_count_fem <- cbind(rep(contig,length(count.fem)),names(count.fem),count.fem)
	#tmp_count_mal <- cbind(rep(contig,length(count.mal)),names(count.mal),count.mal)
	#row.names(tmp_count_fem)<-NULL
	#row.names(tmp_count_mal)<-NULL
	#count.sum.fem <- rbind(count.sum.fem,tmp_count_fem)
	#count.sum.mal <- rbind(count.sum.mal,tmp_count_mal)
}

blasthits <- table(count.blast.hits.mal.norm$V1,count.blast.hits.mal.norm$V4)
rownames(blasthits) <- sub('Backbone_','C_',rownames(blasthits))

#Heat map
my_col <- colorRampPalette(c("green","black","red"))(256)
my_col <- colorRampPalette(c("white","red"))(256)
heatmap.2(blasthits, col= my_col, trace='none', Rowv=FALSE, Colv=FALSE,keysize=1,key.title=NA, key.xlab='Log2FC', cexRow=0.8, cexCol=1, margins=c(5.7,7) , notecol='white', density.info='none')

counts.genes.bp <- read.csv("2018-07-25-tests_FC_normalized_coverage_at_bp_within_genes.txt", 
  sep = "\t", h = T, stringsAsFactors = T)
counts.genes.bp.raw <- counts.genes.bp[which(counts.genes.bp$type == "raw"), ]
counts.genes.bp <- counts.genes.bp[which(counts.genes.bp$type == "norm"), ]

count.blast.hits.norm <- data.frame(count.blast.hits.fem.norm,count.blast.hits.mal.norm$V8,log2((count.blast.hits.fem.norm$V8+1)/(count.blast.hits.mal.norm$V8+1)))
colnames(count.blast.hits.norm) <- c('Contig','start','end','gene','score','sense','query.cover','norm.counts.fem','norm.counts.mal','log2fc')
count.blast.hits.norm$Contig <- sub('Backbone_','C',count.blast.hits.norm$Contig)

write.table(count.blast.hits.norm,'2018-07-27-Blast_female_absent_genes_counts.txt', quote=F, col.names=F,row.names=F, sep='\t')

q <- qplot(gene,Contig,data=count.blast.hits.norm,color=log2fc, size=query.cover, alpha=I(0.5))
q + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_color_gradient2(midpoint=0, low="green", mid="black",
                     high="red", space ="Lab" )

    theme_bw() +
   theme(axis.text.x = element_text(angle = 90, hjust = 1))














