library(seqinr)

file <- read.table('~/Documents/stage_mbelari/results/annotation/2018-06-14-blast_all_genes_vs_hybrid.txt',header=F)

count <- table(file[,1],file[,2])
hpergene <- apply(count,1,sum)

cds <- read.fasta('~/Documents/stage_mbelari/data/raw/Mesorhabditis_belari_JU2817_v2.cds.fa')
cds_names <- getName(cds)
cds_sequences <- getSequence(cds)

genesize <- lengths(cds_sequences)
names(genesize) <- cds_names
file <- data.frame(file,genesize[as.character(file[,1])])
colnames(file)[13]<-'V13'

#genes with one hit
genes_1h <- names(hpergene[hpergene==1])
genes_1h_file <- file[match(genes_1h,file[,1]),]
gene_lenght_covered <- abs(genes_1h_file[,7]-genes_1h_file[,8]) / genes_1h_file[,13]
hist(gene_lenght_covered, breaks=20, main='Query coverage of genes with one hit in blast', xlab='Query coverage percent')

length(which(gene_lenght_covered>0.9))/length(gene_lenght_covered)
#0.8195 genes with one blast hit 
temp <- cbind(genes_1h_file[,3],gene_lenght_covered)
summary(temp[which(gene_lenght_covered>0.9)])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  77.39   99.51  100.00   98.99  100.00  100.00 


