library(seqinr)

genes_list <- c('MBELA.g11678','MBELA.g13063','MBELA.g13483','MBELA.g15528',
'MBELA.g15674','MBELA.g15899','MBELA.g17246','MBELA.g17274','MBELA.g17397',
'MBELA.g17460','MBELA.g20649','MBELA.g21566')

contigs <- c('MBELA.04777','MBELA.06176','MBELA.06579','MBELA.08460','MBELA.08613','MBELA.08772',
'MBELA.10068','MBELA.10105','MBELA.10230','MBELA.10268','MBELA.13121','MBELA.13962')

cds_v2 <- read.fasta('~/Documents/stage_mbelari/data/raw/Mesorhabditis_belari_JU2817_v2.cds.fa')
cds_v2_names <- getName(cds_v2)
cds_v2_seq <- getSequence(cds_v2)

sc_v2 <- read.fasta('~/Documents/stage_mbelari/data/raw/Mesorhabditis_belari_JU2817_v2.scaffolds.fa')
sc_v2_names <- getName(sc_v2)
sc_v2_seq <- getSequence(sc_v2)

pos <- c()

for (i in 1:length(genes_list)){
num <- grep(genes_list[i],cds_v2_names)
pos <- c(pos,num)
}

write.fasta(cds_v2_seq[pos],cds_v2_names[pos],'2018-06-11-13genes_cds_v2.fasta')

posc <- c()

for (i in 1:length(contigs)){
numc <- grep(contigs[i],sc_v2_names)
posc <- c(posc,numc)
}

write.fasta(sc_v2_seq[posc],sc_v2_names[posc],'2018-06-11-13genes_scaffold_v2.fasta')