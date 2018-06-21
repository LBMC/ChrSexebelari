#R functions

getBackbone <- function(list,file){
	library(seqinr)
	genome <- read.fasta('~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta')
	genome_names <- getName(genome)
	genome_seqs <- getSequence(genome)
	list <- paste0('Backbone_',list)
	match_names <- match(list,genome_names)
	#genome_names[match_names]
	#genome_seqs[match_names]
	write.fasta(genome_seqs[match_names],genome_names[match_names],file)
}

getCdsV2 <- function(list,file){
	library(seqinr)
	cds <- read.fasta('~/Documents/stage_mbelari/data/raw/Mesorhabditis_belari_JU2817_v2.cds.fa')
	cds_names <- getName(cds)
	cds_seqs <- getSequence(cds)
	positions <- c()
	for (n in 1:length(list)){
	positions <- c(positions,grep(list[n],cds_names))
	}
	write.fasta(cds_seqs[positions],cds_names[positions],file)
}