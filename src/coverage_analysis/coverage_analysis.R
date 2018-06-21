require(data.table)

contigs.mbela <- fread("~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt", sep = "\t", 
  h = F, stringsAsFactors = F)

cov.female <- fread("~/Documents/stage_mbelari/results/coverage_analysis/2018_06_21_MRDR5_vs_Hybrid_assembly.sorted.filtered.bed", sep = "\t", h = F, stringsAsFactors = F)
cov.male <- fread("~/Documents/stage_mbelari/results/coverage_analysis/2018_06_21_MRDR6_vs_Hybrid_assembly.sorted.filtered.bed", sep = "\t", h = F, stringsAsFactors = F)
# cov.female <- fread("~/Documents/stage_mbelari/results/coverage_analysis/subsample5.bed", sep = "\t", h = F, stringsAsFactors = F)
#cov.male <- fread("~/Documents/stage_mbelari/results/coverage_analysis/subsample6.bed", sep = "\t", h = F, stringsAsFactors = F)
colnames(cov.female) <- c('Contig', 'Reads', 'Length')
colnames(cov.male) <- c('Contig', 'Reads', 'Length')

cov.female <- cov.female[which(cov.female$Contig %in% contigs.mbela$V1), ]
cov.male <- cov.male[which(cov.male$Contig %in% contigs.mbela$V1), ]

source('~/Documents/stage_mbelari/src/coverage_analysis/functions_claire.R')
counts.contigs <- ComputeNormalizedCount(cov.female, cov.male, "contig")

counts.contigs$missing.sexe <- ""; counts.contigs$missing.sexe[which(counts.contigs$counts.raw.female == 0)] <- "female"; 
counts.contigs$missing.sexe[which(counts.contigs$counts.raw.male == 0)] <- "male"; 
write.table(counts.contigs, "results/coverage/FC_normalized_coverage_at_contig.txt", 
    sep = "\t", quote = F, col.names = T, row.names = F)

