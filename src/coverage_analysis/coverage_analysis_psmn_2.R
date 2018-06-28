require(data.table)
library(preprocessCore)

contigs.mbela <- read.table("~/stage_mbelari/results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt", sep = "\t", head = F, row.names=1)

#cov.female.bed <- fread("~/stage_mbelari/results/coverage_analysis/2018_06_21_MRDR5_vs_Hybrid_assembly.sorted.filtered.bed", sep = "\t", h = F, stringsAsFactors = F)
#cov.male.bed <- fread("~/stage_mbelari/results/coverage_analysis/2018_06_21_MRDR6_vs_Hybrid_assembly.sorted.filtered.bed", sep = "\t", h = F, stringsAsFactors = F)

#cov.female.bed <- as.data.frame(cov.female.bed)
#cov.male.bed <- as.data.frame(cov.male.bed)

#reads.f <- tapply(cov.female[,3],cov.female[,1], sum)
#reads.m <- tapply(cov.male[,3],cov.male[,1], sum)

#counts.female <- data.frame(names(reads.f), contigs.mbela[names(reads.f),], reads.f)
#counts.male <- data.frame(names(reads.m), contigs.mbela[names(reads.m),], reads.m)

#colnames(counts.female) <- c('Contig', 'Length', 'Reads')
#colnames(counts.male) <- c('Contig', 'Length', 'Reads')

#write.table(counts.male, '2018_06_21_MRDR6_vs_Hybrid_assembly.sorted.filtered.txt', sep='\t', quote=F, row.names=F, col.names=T)
#write.table(counts.female, '2018_06_21_MRDR5_vs_Hybrid_assembly.sorted.filtered.txt', sep='\t', quote=F, row.names=F, col.names=T)

#Start analysis

cov.female <- fread("~/stage_mbelari/results/coverage_analysis/2018_06_21_MRDR5_vs_Hybrid_assembly.sorted.filtered.txt", sep = "\t", h = T, stringsAsFactors = F)
cov.male <- fread("~/stage_mbelari/results/coverage_analysis/2018_06_21_MRDR6_vs_Hybrid_assembly.sorted.filtered.txt", sep = "\t", h = T, stringsAsFactors = F)

cov.female <- cov.female[which(cov.female$Contig %in% rownames(contigs.mbela)), ]
cov.male <- cov.male[which(cov.male$Contig %in% rownames(contigs.mbela)), ]

source('~/stage_mbelari/src/coverage_analysis/functions_claire.R')
counts.contigs <- ComputeNormalizedCount(cov.female, cov.male, "contig")

counts.contigs$missing.sexe <- ""; counts.contigs$missing.sexe[which(counts.contigs$counts.raw.female == 0)] <- "female";
counts.contigs$missing.sexe[which(counts.contigs$counts.raw.male == 0)] <- "male";
write.table(counts.contigs, "FC_normalized_coverage_at_contig.txt", sep = "\t", quote = F, col.names = T, row.names = F)


