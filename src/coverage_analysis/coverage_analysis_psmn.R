require(data.table)

contigs.mbela <- fread("~/stage_mbelari/results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt", sep = "\t", 
  h = F, stringsAsFactors = F)

cov.female <- fread("~/stage_mbelari/results/coverage_analysis/2018_06_21_MRDR5_vs_Hybrid_assembly.sorted.filtered.bed", sep = "\t", h = F, stringsAsFactors = F)
cov.male <- fread("~/stage_mbelari/results/coverage_analysis/2018_06_21_MRDR6_vs_Hybrid_assembly.sorted.filtered.bed", sep = "\t", h = F, stringsAsFactors = F)


counts.male <- contigs.mbela
counts.male <- data.frame(counts.male,vector(length=nrow(counts.male)))

for (n in 1:nrow(contigs.mbela)){
con <- contigs.mbela[n,1]
counts.male[n,3] <- sum(cov.male[which(cov.male[,1]==con),3])
}

counts.female <- contigs.mbela
counts.female <- data.frame(counts.female,vector(length=nrow(counts.female)))

for (n in 1:nrow(contigs.mbela)){
con <- contigs.mbela[n,1]
counts.female[n,3] <- sum(cov.female[which(cov.female[,1]==con),3])
}

cov.female <- counts.female
cov.male <- counts.male

colnames(cov.female) <- c('Contig', 'Length', 'Reads')
colnames(cov.male) <- c('Contig', 'Length', 'Reads')

cov.female <- cov.female[which(cov.female$Contig %in% contigs.mbela$V1), ]
cov.male <- cov.male[which(cov.male$Contig %in% contigs.mbela$V1), ]

source('~/stage_mbelari/src/coverage_analysis/functions_claire.R')
counts.contigs <- ComputeNormalizedCount(cov.female, cov.male, "contig")

counts.contigs$missing.sexe <- ""; counts.contigs$missing.sexe[which(counts.contigs$counts.raw.female == 0)] <- "female"; 
counts.contigs$missing.sexe[which(counts.contigs$counts.raw.male == 0)] <- "male"; 
write.table(counts.contigs, "FC_normalized_coverage_at_contig.txt", sep = "\t", quote = F, col.names = T, row.names = F)

