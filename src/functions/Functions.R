#R functions
#to call source('~/Documents/stage_mbelari/src/functions/Functions.R')

#To install: install.packages('seqinr')

#Creates a fasta file with selected Backbone sequences 
#list is the list of the number of Backbones:  c(11,495,601)
#file is the name of the fasta file (add .fasta)
#getBackbone(c(11,495),'name_file.fasta')
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

#Creates a fasta file with selected Backbone sequences from selected regions (not the complete Backbone)
#list is the list of the number of Backbones:  c(11,495,601)
#file is the name of the fasta file (add .fasta)
#start is the list of start positions for each backbone
#end is the list of end positions for each backbone
frag.Backbone <- function(list,file,start,end){
	library(seqinr)
	genome <- read.fasta('~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta')
	genome_names <- getName(genome)
	genome_seqs <- getSequence(genome)
	list <- paste0('Backbone_',list)
	match_names <- match(list,genome_names)
	tmp <- genome_seqs[match_names]
	for (n in 1:length(list)){
	tmp[[n]] <- tmp[[n]][start[n]:end[n]]	
	}
	write.fasta(tmp,genome_names[match_names],file)
}

#Creates a fasta file with Cds sequences from illumina genome. 
#list is the list of the genes: c(g110,g1542)
#file is the name of the fasta file (add .fasta)
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

#Creates a fasta file with Cds sequences from hybrid assembly.
#list is the list of the genes c(g110,g1542)
#file is the name of the fasta file (add .fasta)
getCdsH <- function(list,file){
	library(seqinr)
	cds <- read.fasta('~/Documents/stage_mbelari/results/annotation/BRAKER/soft_masked/augustus.hints.codingseq')
	cds_names <- getName(cds)
	cds_seqs <- getSequence(cds)
	tmp_names <- sub('Backbone_[0-9]+\\.','',getName(cds))
	tmp_names <- sub('\\.t[0-9]+','',tmp_names)
	match_names <- match(list,tmp_names)
	write.fasta(cds_seqs[match_names],cds_names[match_names],file)
}
