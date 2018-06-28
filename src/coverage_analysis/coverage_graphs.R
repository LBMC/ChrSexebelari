require(data.table)
source('http://bioconductor.org/biocLite.R')
biocLite('preprocessCore')
library(preprocessCore)

source('~/Documents/stage_mbelari/src/coverage_analysis/functions_claire.R')

contigs.mbela <- read.table("~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt", sep = "\t", head = F, row.names=1)
cov.female <- fread("~/Documents/stage_mbelari/results/coverage_analysis/2018_06_21_MRDR5_vs_Hybrid_assembly.sorted.filtered.bed", sep = "\t", h = T, stringsAsFactors = F)
cov.male <- fread("~/Documents/stage_mbelari/results/coverage_analysis/2018_06_21_MRDR6_vs_Hybrid_assembly.sorted.filtered.bed", sep = "\t", h = T, stringsAsFactors = F)

cov.female[cov.female[,1]==contigs.mbela[1,1],]


colnames(cov.female) <- c('Contig','Length','Reads','Unmapped')
colnames(cov.male) <- c('Contig','Length','Reads','Unmapped')

cov.female <- cov.female[which(cov.female$Contig %in% rownames(contigs.mbela)), ]
cov.male <- cov.male[which(cov.male$Contig %in% rownames(contigs.mbela)), ]

counts.contigs <- ComputeNormalizedCount(cov.female, cov.male, "contig")

counts.contigs$missing.sexe <- ""; counts.contigs$missing.sexe[which(counts.contigs$counts.raw.female == 0)] <- "female";
counts.contigs$missing.sexe[which(counts.contigs$counts.raw.male == 0)] <- "male";
write.table(counts.contigs, "2018-06-21-FC_normalized_coverage_at_contig.txt", sep = "\t", quote = F, col.names = T, row.names = F)


