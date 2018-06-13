library(seqinr)

#read files

#blast output
file <- read.table('~/Documents/stage_mbelari/results/13genes/2018-06-01-13genes_in_hybrid_genome.txt',header=F)
sense <- file[,9]>file[,10]
file <- data.frame(file,sense)

#genome contigs
genome <- read.fasta('~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta')
genome_names <- getName(genome)
genome_seqs <- getSequence(genome)

#13 genes
genes_13 <- read.fasta('~/Documents/stage_mbelari/data/raw/2018-05-18-Complete_13_genes.fasta')
genes_13_names <- getName(genes_13)
genes_sequences <- getSequence(genes_13)

#cds seqs
cds_13 <- read.fasta('~/Documents/stage_mbelari/data/raw/mbela_cds_JU2817.fa')
cds_names <- getName(cds_13)
cds_sequences <- getSequence(cds_13)

genes_tested <- unique(names(table(file[,1])))

outer_regions <- function(gene,file,nb,th){
	subfile <- file[file[,1]==gene,]
	scores <- subfile[,12]/subfile[1,12]
	subfile <- data.frame(subfile,scores)
	subfile <- subfile[subfile[,14]>th,]
	reads <- subfile[,2]
	sizes <- lengths(genome_seqs[match(reads,genome_names)])

	bgns <- c()
	ends <- c()
	for (x in 1:length(reads)){
		read <- as.character(reads[x])
		if (subfile[x,9]>nb){
			bgns[x] <- subfile[x,9]-nb
		} else {
			bgns[x] <- 1
		}
		if ((sizes[x]-subfile[x,10])>nb){
			ends[x] <- subfile[x,10]+nb
		} else {
			ends[x] <- sizes[x]
		}

		}
		
	seqs_to_cut <- genome_seqs[match(reads,genome_names)]
	finalseqs <- list()
	for (s in 1:length(seqs_to_cut)){
		finalseqs[[s]] <- seqs_to_cut[[s]][bgns[s]:ends[s]]
		if (sense[s]==TRUE){
			 finalseqs[[s]] <- rev(comp(finalseqs[[s]]))
		}
	}
	
	finalseqs[[s+1]] <- genes_sequences[genes_13_names==gene][[1]]
	finalseqs[[s+2]] <- cds_sequences[[grep(gene,cds_names)]]
	finalnames <- genome_names[match(reads,genome_names)]
	finalnames <- c(finalnames, gene, paste0(gene,'_cds'))
	
	namefile <- paste0('Hybrid_contigs_w_',gene,'.txt')
	write.fasta(finalseqs,finalnames,file.out=namefile)	
}

for (l in (1:length(genes_tested))){
	gene <- genes_tested[l]
	outer_regions(gene,file,0,0.75)
	}

count_genes_in_scaffolds <- table(file[,1],file[,2])
#write.csv(count_genes_in_scaffolds,'Genes_in_scaffolds.csv')

Large_contig <- file[file[,2]=='Backbone_495',]
Check <- c(1,1,1,0,1,0,1)
Large_contig <- data.frame(Large_contig,Check)
Large_contig[Check==1,]

genes_13 <- read.fasta('~/Documents/stage_mbelari/data/raw/2018-05-18-Complete_13_genes.fasta')
cds_names <- sub('\\.t[0-9]*','',cds_names)
getcds <- match(getName(genes_13)[-1],cds_names)

file <- read.table('~/Documents/stage_mbelari/results/13genes/2018-06-01-13genes_in_hybrid_genome.txt',header=F)
sense <- file[,9]>file[,10]
file <- data.frame(file,sense)
count_genes_in_scaffolds <- table(file[,1],file[,2])
for (l in 1:length(table(file[,1]))){
	subfile <- file[file[,1]==rownames(count_genes_in_scaffolds)[l],]
	temp <- subfile[,12]
	count_genes_in_scaffolds[subfile[1,1],subfile[,2]] <- subfile[,12]
	repeated <- which(table(names(temp))>1)
	if(length(repeated)!=0){
	count_genes_in_scaffolds[subfile[1,1],names(repeated)] <- max(subfile[subfile[2]==names(repeated),12])
	}
}
final_count <- apply(count_genes_in_scaffolds,1,function(x){x/max(x)})
library(RColorBrewer)
#my_col=colorRampPalette(brewer.pal(8, "Blues"))(10)
my_col=colorRampPalette(c("white","purple"))(10)
heatmap(final_count, col= my_col)
library(gplots)
heatmap.2(final_count, col= my_col, trace='none', Rowv=FALSE, Colv=FALSE)
png('Genes_in_contigs.png', width=1000,height=1000)
heatmap.2(final_count, col= my_col, trace='none', Rowv=FALSE, Colv=FALSE,keysize=1,key.title=NA, key.xlab='Alignment score', cexRow=0.9, cexCol=1.3, margins=c(5.5,6.8))
dev.off()

library(genoPlotR)

dna_seg_list <- list()

#drawcontigs <- names(which(table(file[,2])>1))
drawcontigs <- c('Backbone_495','Backbone_1072','Backbone_183')
#drawcontigs <- c('Backbone_1052','Backbone_1092','Backbone_38')
annot_list <- list()
for (l in 1:length(drawcontigs)){
	subfilea <- file[file[,2]==drawcontigs[l],]
	gene_names <- as.character(subfilea[,1])
	starts <- subfilea[,9]
	ends <- subfilea[,10]
	strands <- as.numeric(subfilea[,13])
	strands[strands==0] <- -1
	cols <- colorRampPalette(c("darksalmon","firebrick4"))(length(gene_names))
	df <- data.frame(name=gene_names,start=starts,end=ends,strand=strands,col=cols)
	dna_seg_list[[l]] <- as.dna_seg(df, fill=cols, lty=2)
	annot_list[[l]] <- annotation(middle(dna_seg_list[[l]]), text=gene_names, rot=0, col=dna_seg_list[[l]]$col)
}

col <- c('brown3','deepskyblue3','forestgreen','darkslateblue','coral2','darkgoldenrod1','darkolivegreen2','antiquewhite3','cornflowerblue')
col1 <- col[c(1,1,2,2,3,3,4)]
col2 <- col[c(5,1,2,4)]
col3 <- col[c(6,7,8)]
dna_seg_list[[1]]$col <-col1
dna_seg_list[[2]]$col <-col2
dna_seg_list[[3]]$col <-col3
dna_seg_list[[1]]$fill <-col1
dna_seg_list[[2]]$fill <-col2
dna_seg_list[[3]]$fill <-col3
png('Contigs.png', width=800, height=300)
#png('Contigs.png')
#pdf('Contigs.pdf', width=1000, height=1000)
plot_gene_map(dna_seg_list, dna_seg_labels=drawcontigs, dna_seg_labes_cex=1.2,annotations=annot_list, annotations.cex=1.5)
dev.off()

