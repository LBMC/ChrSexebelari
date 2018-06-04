library(seqinr)

#read files

#blast output
file <- read.table('~/Documents/stage_mbelari/results/13genes/2018-06-01-13genes_in_hybrid_genome.txt',header=F)
sense <- file[,9]>file[,10]
file <- data.frame(file,sense)

#genome contigs
genome <- read.fasta('~/Documents/stage_mbelari/results/hybrid_test/2018-06-01-DGB2OLC/final_assembly.fasta')
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

outer_regions <- function(gene,file){
	subfile <- file[file[,1]==gene,]
	reads <- subfile[,2]
	sizes <- lengths(genome_seqs[match(reads,genome_names)])

	bgns <- c()
	ends <- c()
	for (x in 1:length(reads)){
		read <- as.character(reads[x])
		if (subfile[x,9]>500){
			bgns[x] <- subfile[x,9]-500
		} else {
			bgns[x] <- 1
		}
		if ((sizes[x]-subfile[x,10])>500){
			ends[x] <- subfile[x,10]+500
		} else {
			ends[x] <- sizes[x]
		}

		}
		
	seqs_to_cut <- genome_seqs[match(reads,genome_names)]
	for (s in 1:length(seqs_to_cut)){
		finalseqs[[s]] <- seqs_to_cut[[s]][bgns[s]:ends[s]]
		if (sense[s]==TRUE){
			 finalseqs[[s]] <- rev(comp(finalseqs[[s]]))
		}
	}
	namefile <- paste0('Contigs_w_',gene,'.txt')
	write.fasta(finalseqs,genome_names[match(reads,genome_names)],file.out=namefile)	
}

for (l in (1:genes_tested)){
	gene <- genes_tested[l]
	outer_regions(gene,file)
	}