library(seqinr)

getBackbone <- function(list,file){
	library(seqinr)
	genome <- read.fasta('~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta')
	genome_names <- getName(genome)
	genome_seqs <- getSequence(genome)
	list <- paste0('Backbone_',list)
	match_names <- match(list,genome_names)
	genome_names[match_names]
	genome_seqs[match_names]
	write.fasta(genome_seqs[match_names],genome_names[match_names],file)
}

#genome contigs
genome <- read.fasta('~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta')
genome <- read.fasta('~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-07/final_assembly_min500_001.fasta')
genome_names <- getName(genome)
genome_seqs <- getSequence(genome)

table <- read.csv('Genes_in_scaffolds.csv')
male_backbones <- names(table)[which(table[13,]=='m')]

males_names <- match(male_backbones,genome_names)
genome_names[males_names]
genome_seqs[males_names]

write.fasta(genome_seqs[males_names],genome_names[males_names],'2018-06-08-Contigs_males.fasta')

L50<-length(genome_seqs)/2
get <- order(lengths(genome_seqs))<=L50
sub50_seqs <- genome_seqs[get]
sub50_names <- genome_names[get]
write.fasta(sub50_seqs,sub50_names,'2018-06-12-sub_50_hybrid_3_500.fasta')

