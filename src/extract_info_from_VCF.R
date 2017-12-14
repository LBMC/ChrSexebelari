require(data.table)
require(dplyr)
library(foreach)
library(doParallel)
library(gridExtra)

# size of contigs
size.genome.mbelari <- read.table("data/ReferenceGenomes/2017_09_13_Mbelari.sizes.genome", h = F, sep = "\t", stringsAsFactors = F)
# annotation
gff.mbelari <- read.csv("data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2.gff3", h = F, sep = "\t", stringsAsFactors = F)
# genes annotation
gff.mbelari.genes <- gff.mbelari[which(gff.mbelari$V3 == "gene"), ]

do.computation.per.pool <- F
do.summary.raw <- F

##### Extract informations from VCF
if (do.computation.per.pool) {
  # Unzip vcf files
  system("gunzip -c results/call_var/2017_09_21_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz > results/call_var/2017_09_21_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf")
  system("gunzip -c results/call_var/2017_09_21_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz > results/call_var/2017_09_21_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf")
  
  variants.female <- fread("results/call_var/2017_09_21_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf", sep = "\t",  sep2 = ";", h = T, skip = 19286)
  variants.male <- fread("results/call_var/2017_09_21_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf", sep = "\t",  sep2 = ";", h = T, skip = 19286)

  # Get nb of alleles per pool (AN), number of reads "equal" to ref/number of reads equal to alt at a given position (DP4 field)
  info.female <- ExtractCountFromVCF(variants.female, "female")
  info.male <- ExtractCountFromVCF(variants.male, "male")
  
  variants.female <- cbind(variants.female, info.female)
  variants.male  <- cbind(variants.male, info.male)
  colnames(variants.female)[1] <- "CHROM"
  colnames(variants.male)[1] <- "CHROM"
  
  write.table(variants.female, "results/call_var/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts.txt", sep = "\t",  quote = F, row.names = F)
  system("bash src/date.sh MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts.txt")
  write.table(variants.male, "results/call_var/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts.txt", sep = "\t",  quote = F, row.names = F)
  system("bash src/date.sh results/call_var/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts.txt")
  system("bash src/date.sh results/call_var/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts.txt")
}

tab.variants.female <- fread("results/call_var/2017_09_25_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts.txt", sep = "\t", h = T, stringsAsFactors = F)
tab.variants.male <- fread("results/call_var/2017_09_25_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts.txt", sep = "\t", h = T, stringsAsFactors = F)

##### Summary of detected variants
if(do.summary.raw) {
	summary.raw.variants <- NULL
	for (sexe in c("male", "female")) {
		eval(parse(text = paste("tmp = tab.variants.", sexe, sep = "")))
		tmp.tab <- table(tmp$is.INDEL)
		nb.SNPs.all <- tmp.tab["0"]
		nb.INDELs.all <- tmp.tab["1"]
		nb.pluriallelic.sites <- length(which(tmp$nb.alleles>2))
		eval(parse(text = paste("nb.median.depth.all <- median(tmp$tot.", sexe, ")", sep = "")))
		eval(parse(text = paste("nb.mean.depth.all <- mean(tmp$tot.", sexe, ")", sep = "")))
		eval(parse(text = paste("nb.median.freq.all <- median(100*tmp$count.alt.", sexe, "/tmp$tot.", sexe, ")", sep = "")))
		eval(parse(text = paste("nb.mean.freq.all <- mean(100*tmp$count.alt.", sexe, "/tmp$tot.", sexe, ")", sep = "")))
	
		##### Compute what is in genic areas only
		cores=detectCores()
		cl <- makeCluster(cores[1]-1) 
		registerDoParallel(cl)
		all <- foreach(ctg=unique(tmp$CHROM), .combine=rbind) %dopar% {
			temp <- tmp[which(tmp$CHROM == ctg), ]
			temp.gff <- gff.mbelari.genes[which(gff.mbelari.genes$V1 == ctg), ]
			final <- cbind(temp, data.frame(genes = sapply(temp$POS, function(x) paste(temp.gff$V9[which(temp.gff$V4<= x & temp.gff$V5>=x)], collapse = "_")), stringsAsFactors = F))
			final
		}
		stopCluster(cl)
		write.table(all, paste("results/call_var/", sexe, "_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts_genic_information.txt", sep = ""), sep = "\t",  quote = F, row.names = F)
		system(paste("bash src/date.sh results/call_var/", sexe, "_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts_genic_information.txt", sep = ""))

		all.genic <- all[-which(all$genes == ""), ]
		tmp.tab<- table(all.genic$is.INDEL)
		nb.SNPs.genic <- tmp.tab["0"]
		nb.INDELs.genic <- tmp.tab["1"]
		eval(parse(text = paste("nb.median.depth.genic <- median(all.genic$tot.", sexe, ")", sep = "")))
		eval(parse(text = paste("nb.mean.depth.genic <- mean(all.genic$tot.", sexe, ")", sep = "")))
		eval(parse(text = paste("nb.median.freq.genic <- median(100*all.genic$count.alt.", sexe, "/all.genic$tot.", sexe, ")", sep = "")))
		eval(parse(text = paste("nb.mean.freq.genic <- mean(100*all.genic$count.alt.", sexe, "/all.genic$tot.", sexe, ")", sep = "")))

		##### Compute distance to the nearest INDELs
		cores=detectCores()
		cl <- makeCluster(cores[1]-1) #not to overload your computer
		registerDoParallel(cl)
		min.dist <- foreach(i=1:dim(all.genic)[1], .combine=c) %dopar% {
			tmp <- all.genic[i, ]
			if(tmp$is.INDEL == 1) {
				temp <- all.genic[which(all.genic$CHROM == tmp$CHROM & all.genic$POS != tmp$POS), ]
				d <- temp$POS-tmp$POS
				dist <- d[which(abs(d) == min(abs(d)))]
			} else{dist <- NA}
			dist
		}
		stopCluster(cl)
	
		##### Distance from INDEL to nearest variant
		pdf(paste("results/call_var/distance_from_INDEL_to_nearest_variant_", sexe, ".pdf", sep = ""))
		hist(na.omit(min.dist), breaks = 40000, xlim = c(-20,20), main = paste("distance_from_INDEL_to_nearest_variant ", sexe, sep = ""), xlab = "distance in bp")
		dev.off()	
		system(paste("bash src/date.sh results/call_var/distance_from_INDEL_to_nearest_variant_", sexe, ".pdf", sep = ""))

		summary.raw.variants <- rbind(summary.raw.variants, data.frame(sexe = sexe, nb.SNPs.all = nb.SNPs.all, nb.INDELs.all = nb.INDELs.all, nb.pluriallelic.sites = nb.pluriallelic.sites, nb.median.depth.all = nb.median.depth.all, nb.mean.depth.all = nb.mean.depth.all, median.freq.all = nb.median.freq.all, mean.freq.all = nb.mean.freq.all, nb.SNPs.genic = nb.SNPs.genic, nb.INDELs.genic = nb.INDELs.genic, nb.median.depth.genic = nb.median.depth.genic, nb.mean.depth.genic = nb.mean.depth.genic, median.freq.genic = nb.median.freq.genic, mean.freq.genic = nb.mean.freq.genic, stringsAsFactors = F))	
	}
	write.table(summary.raw.variants, "results/call_var/summary_raw_variants.txt", sep = "\t",  quote = F, row.names = F)
	system("bash src/date.sh results/call_var/summary_raw_variants.txt")
}

s <- read.table("results/call_var/summary_raw_variants.txt", sep = "\t", h = T)
pdf("results/call_var/summary_raw_variants.pdf", height=4, width=25)
grid.table(s)
dev.off()
system("bash src/date.sh results/call_var/summary_raw_variants.pdf")


##### Summary of density of SNPs for raw variants per contig 
nb.test.female <- sapply(size.genome.mbelari$V1, function(x) length(which(tab.variants.female$CHROM == x & tab.variants.female$is.INDEL == 0)))
nb.test.male <- sapply(size.genome.mbelari$V1, function(x) length(which(tab.variants.male$CHROM == x & tab.variants.male$is.INDEL == 0)))
density.tests.sexe <- data.frame(contig = size.genome.mbelari$V1, length = size.genome.mbelari$V2, nb.test.female= nb.test.female, nb.test.male = nb.test.male, stringsAsFactors = F)
density.tests.sexe$den.male <- density.tests.sexe$nb.test.male/density.tests.sexe$length
density.tests.sexe$den.female <- density.tests.sexe$nb.test.female/density.tests.sexe$length
log2.den.var <- log2(density.tests.sexe$den.female/density.tests.sexe$den.male)
pdf("results/call_var/log2_density_ratio_raw_SNPs.pdf")
hist(log2.den.var, main = "log2(density SNP female/density SNP male) on raw detected SNPs", xlab  = "log2(density SNP female/density SNP male)", breaks = 100)
abline(v = 0, col = "black", lty = 2)
dev.off()
system("bash src/date.sh results/call_var/log2_density_ratio_raw_SNPs.pdf")

