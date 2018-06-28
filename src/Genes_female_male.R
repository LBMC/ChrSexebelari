library(seqinr)

source('~/Documents/stage_mbelari/src/Functions.R')

FCgenes <- read.table('2017_12_06_FC_normalized_coverage_at_gene.txt', head=T, sep='\t')
male_absent <- FCgenes[FCgenes$missing.sexe=='male',]
genes_male_absent <- sub('MBELA\\.','',as.character(male_absent[,1]))

#getCdsV2(genes_male_absent,'List_genes_absent_in_males.txt')
female_absent <- FCgenes[FCgenes$missing.sexe=='female',]

top_50_fem_abs <- FCgenes[order(FCgenes$log2.norm.FC)[1:50],]
genes_top_50_fem_abs <- sub('MBELA\\.','',as.character(top_50_fem_abs[,1]))

getCdsV2(genes_top_50_fem_abs,'2018-06-20-Top_50_genes_fem_abs.fasta')

Blast_results <- read.table('~/Documents/stage_mbelari/results/13genes/version2/2018-06-20-Genes_absent_females.txt',head=FALSE)
genes_in_backbones <- table(Blast_results[,1],Blast_results[,2])

genes_in_backbones[genes_in_backbones>=1]<-1
apply(genes_in_backbones,2,sum)

