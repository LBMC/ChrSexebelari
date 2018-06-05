library(seqinr)

#read files

#blast output
file <- read.table('~/Documents/stage_mbelari/results/13genes/2018-06-05-13genes_in_illumina_genome.txt',header=F)
sense <- file[,9]>file[,10]
file <- data.frame(file,sense)

#genome contigs
genome <- read.fasta('~/Documents/stage_mbelari/data/raw/mbela_assembly_JU2817.fa')
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
	
	namefile <- paste0('Illumina_contigs_w_',gene,'.txt')
	write.fasta(finalseqs,finalnames,file.out=namefile)	
}

for (l in (1:length(genes_tested))){
	gene <- genes_tested[l]
	outer_regions(gene,file,0,0.75)
	}

count_genes_in_scaffolds <- table(file[,1],file[,2])
write.csv(count_genes_in_scaffolds,'Genes_in_scaffolds.csv')

Large_contig <- file[file[,2]=='Backbone_495',]
Check <- c(1,1,1,0,1,0,1)
Large_contig <- data.frame(Large_contig,Check)
Large_contig[Check==1,]

