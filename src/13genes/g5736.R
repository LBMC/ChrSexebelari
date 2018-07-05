library(seqinr)
source('~/Documents/stage_mbelari/src/Functions.R')
#read files

#blast output
file <- read.table('~/Documents/stage_mbelari/results/13genes/version2/2018-06-28-Genes_absent_females.txt',header=F)

#genome contigs
genome <- read.fasta('~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta')
genome_names <- getName(genome)
genome_seqs <- getSequence(genome)

#13 genes
#MBELA.g5736.t1

subfile <- file[file[,1]=='MBELA.g5736.t1',]
backbones <- unique(as.character(subfile[,2]))
getBackbone(sub('Backbone_','',backbones), 'Backbones.g5736.t1')




