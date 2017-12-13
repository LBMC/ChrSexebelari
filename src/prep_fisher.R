require(data.table)
require(dplyr)
library(foreach)
library(doParallel)
library(gridExtra)

##### Select contigs with size >= 1,000 bp
size.genome.mbelari <- read.table("data/ReferenceGenomes/2017_09_13_Mbelari.sizes.genome", h = F, sep = "\t", stringsAsFactors = F)
size.contigs1000 <- size.genome.mbelari$V1[which(size.genome.mbelari$V2 >= 1000)]

do.prep.fisher <- F

#### Pick information in variants
tab.variants.female <- fread("results/call_var/2017_11_09_female_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts_genic_information.txt", sep = "\t", h = T, stringsAsFactors = F)
tab.variants.male <- fread("results/call_var/2017_11_09_male_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts_genic_information.txt", sep = "\t", h = T, stringsAsFactors = F)

if (do.prep.fisher) {
 
  ##### Variants with same position in both pools
  merge.sexe.common <- merge(tab.variants.male[, c("CHROM", "POS", "REF", "ALT", "QUAL", "nb.alleles", "count.ref.male", "count.alt.male", "tot.male", "is.INDEL", "genes")], tab.variants.female[, c("CHROM", "POS", "REF", "ALT", "QUAL", "nb.alleles", "count.ref.female", "count.alt.female", "tot.female", "is.INDEL", "genes")], by = c("CHROM",  "POS", "REF", "genes"))

  ##### All variants in both pools 
  merge.sexe.all <- merge(tab.variants.male[, c("CHROM", "POS", "REF", "ALT", "QUAL", "nb.alleles", "count.ref.male", "count.alt.male", "tot.male", "is.INDEL","genes")], tab.variants.female[, c("CHROM", "POS", "REF", "ALT", "QUAL", "nb.alleles", "count.ref.female", "count.alt.female", "tot.female", "is.INDEL", "genes")], by = c("CHROM",  "POS", "REF", "genes"), all.x = T, all.y = T)

  ##### Variants with same position and same alt allele in both pools
  merge.sexe.common.biallelic <- merge(tab.variants.male[, c("CHROM", "POS", "REF", "ALT", "QUAL", "nb.alleles", "count.ref.male", "count.alt.male", "tot.male","is.INDEL", "genes")], tab.variants.female[, c("CHROM", "POS", "REF", "ALT", "QUAL", "nb.alleles", "count.ref.female", "count.alt.female", "tot.female", "is.INDEL", "genes")], by = c("CHROM",  "POS", "REF", "ALT", "genes"))
 
  ##### Apply filter on contig size, DP (>1 in total count) and distance to INDEL:
  # for variants with same position and same alt allele in both pools
  tmp.common.biallelic <- merge.sexe.common.biallelic[which(merge.sexe.common.biallelic$tot.male> 1 | merge.sexe.common.biallelic$tot.female>1), ]
  tmp.common.biallelic$ID <- 1:dim(tmp.common.biallelic)[1]

  # for all variants in both pools
  tmp.all <- merge.sexe.all[which(merge.sexe.all$tot.male>1 | merge.sexe.all$tot.female>1), ]
  tmp.all$ID <- 1:dim(tmp.all)[1]

  ##### Compute distance to the nearest INDELs to remove SNP at 1bp away from INDEL and remove INDELs
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  remove.dist <- foreach(i=1:dim(tmp.all)[1], .combine=c) %dopar% {
  	tmp2 <- tmp.all[i, ]
  	if(tmp2$is.INDEL == 0) {
		temp <- filter(tmp.all, CHROM == tmp2$CHROM & tmp.all$is.INDEL == 1 & tmp.all$ID != tmp2$ID)
		if(dim(temp)[1]>0){
			d <- temp$POS-tmp2$POS
			dist <- min(abs(d))
			if (dist <= 1){
				remove <- tmp2$ID
			}else{
				remove <- NA
			}
		}else{
			remove <- NA
		}
	} else{
		remove <- NA
	}
	remove
  }
  stopCluster(cl)

  ##### Remove SNP at 1bp away from INDEL
  ids.remove <- remove.dist[-which(is.na(remove.dist))]
  merge.sexe.all.filt <- tmp.all[-which(tmp.all$ID %in% ids.remove), ]
  summ3b <- SummarizeSNPsINDELsWithinDataFrame(merge.sexe.all.filt, suff = "all position, DP>1, dist.to.INDEL>1") 
  summary.filt <- rbind(summ1, summ2, summ3, summ1b, summ2b, summ3b)

  ##### Add same.ALT column if ALT allele differs between male and female
  merge.sexe.all.filt$same.ALT <- 1
  merge.sexe.all.filt$same.ALT[which((merge.sexe.all.filt$ALT.x == merge.sexe.all.filt$ALT.y) == F)] <- 0

  ##### Add contig length
  mat.size.ctg <- matrix(size.genome.mbelari$V2, nrow = 1, dimnames = list(NULL, size.genome.mbelari$V1))
  merge.sexe.all.filt$contig_length <- mat.size.ctg[,  merge.sexe.all.filt$CHROM]

  ##### Summary of filtered variants
  tab.snp <- table(merge.sexe.all.filt$same.ALT[which(merge.sexe.all.filt$is.INDEL == 0)])
  tab.indel <- table(merge.sexe.all.filt$same.ALT[which(merge.sexe.all.filt$is.INDEL == 1)])
  pdf("results/call_var/summary_alleles_ALT_between_pools.pdf", height=4, width=10)
  par(mfrow = c(1, 2))
  bb.snp <- barplot(tab.snp, names.arg = c("different ALT allele", "same ALT allele"), main = "SNP identity between both pools", cex.arg = 0.75)
  text(bb.snp, max(tab.snp)/2,labels=tab.snp) 
  bb.indel <- barplot(tab.indel, names.arg = c("different ALT allele", "same ALT allele"), main = "INDELs identity between both pools", cex.arg = 0.75) 
  text(bb.indel, max(tab.indel)/2,labels=tab.indel) 
  dev.off()
  system("bash src/date.sh results/call_var/summary_alleles_ALT_between_pools.pdf")

  write.table(summary.filt, "results/call_var/summary_filt_variants.txt", sep = "\t",  quote = F, row.names = F)
  write.table(merge.sexe.all.filt, "results/call_var/merge.sexe.all.filt.txt", sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(merge.sexe.all.filt[which(merge.sexe.all.filt$is.INDEL == 0),], "results/call_var/merge.sexe.all.filt.SNP.txt", sep = "\t", row.names = F, col.names = T, quote = F)

  system("bash src/date.sh results/call_var/summary_filt_variants.txt")
  system("bash src/date.sh results/call_var/merge.sexe.all.filt.txt")
  system("bash src/date.sh results/call_var/merge.sexe.all.filt.SNP.txt")
}

f <- read.csv("results/call_var/2017_12_05_summary_filt_variants.txt", sep = "\t", h = T)
pdf("results/call_var/summary_filt_variants.pdf", height=4, width=12)
grid.table(f)
dev.off()
system("bash src/date.sh results/call_var/summary_filt_variants.pdf")

##### Summary of density of SNPs for filtered variants per contig 
tests  <- fread("results/call_var/2017_11_13_merge.sexe.all.filt.SNP.txt", sep = "\t", h = T)
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



# Plot histogram of frequency of filtered variants 
data.filt.SNPs <- read.csv("results/call_var/2017_11_09_merge.sexe.all.filt.SNP.txt", sep = "\t", h = T, stringsAsFactors = F)

pdf("results/call_var/2017_11_13_merge.sexe.common.biallelic.filt.SNP.frequency.male.female.pool.pdf")
hist(data.filt.SNPs$count.alt.male/data.filt.SNPs$tot.male, col = adjustcolor("blue", 0.75), breaks = 20, xlim = c(0,1), xlab = "Variants frequency", main = "Histogram of filtered SNPs frequencies")
hist(data.filt.SNPs$count.alt.female/data.filt.SNPs$tot.female, col = adjustcolor("red", 0.75), add = T, breaks = 20, xlim = c(0,1))
legend("topleft", c("male variants frequency", "female variants frequency"), fill = c("blue", "red"))
dev.off()
