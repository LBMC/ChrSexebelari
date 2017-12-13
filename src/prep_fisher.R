require(data.table)
require(dplyr)
library(foreach)
library(doParallel)
library(gridExtra)

##### Select contigs with size >= 1,000 bp
size.genome.mbelari <- read.table("data/ReferenceGenomes/2017_09_13_Mbelari.sizes.genome", h = F, sep = "\t", stringsAsFactors = F)

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

  # Manage is.INDEL when variants is found only once (NA) or found in both but with different type (FALSE)
  test <- apply(tmp.all[, paste("is.INDEL.", c("x", "y"), sep = "")], 1, function(x) ifelse(is.na(x[["is.INDEL.x"]]) | is.na(x[["is.INDEL.y"]]), NA, 
x[["is.INDEL.x"]] == x[["is.INDEL.y"]]))

  tab.test <- table(test, useNA = "always") 
  tmp.all$same_class_variant <- test

  tmp.all.SNP <- tmp.all[which(tmp.all$is.INDEL.x == 0 & tmp.all$is.INDEL.y == 0), ]
  tmp.all.INDEL <- tmp.all[which(tmp.all$is.INDEL.x == 1 & tmp.all$is.INDEL.y == 1), ] 
  tmp.all.SNP.in.one.sexe <- tmp.all[which(tmp.all$same_class_variant == F), ] 
  tmp.all.in.one.sexe <- tmp.all[which(is.na(tmp.all$same_class_variant)), ] 

  # Add supplementary SNP for variants that are both SNP and INDEL
  add.SNP.in.one.sexe <- NULL
  for(i in 1:dim(tmp.all.SNP.in.one.sexe)[1]){
	temp <- as.data.frame(tmp.all.SNP.in.one.sexe[i, ])
        sex.indel <- c("x", "y")[which(temp[, c("is.INDEL.x", "is.INDEL.y")] == 1)]
	sex <- ifelse(sex.indel == "x", "male", "female")
	change <- paste(c("is.INDEL.", "ALT.", "QUAL.", "is.INDEL.", "nb.alleles."), sex.indel, sep = "")
	for (c in change) {temp[1, c]<- NA}
	change <- paste(c("count.ref.", "count.alt.", "tot."), sex, sep = "")
	for (c in change) {temp[1, c]<- 0}
	add.SNP.in.one.sexe <- rbind(add.SNP.in.one.sexe, temp)
  }

  # Add supplementary INDEL for variants that are both SNP and INDEL
  add.INDEL.in.one.sexe <- NULL
  for(i in 1:dim(tmp.all.SNP.in.one.sexe)[1]){
	temp <- as.data.frame(tmp.all.SNP.in.one.sexe[i, ])
        sex.indel <- c("x", "y")[which(temp[, c("is.INDEL.x", "is.INDEL.y")] == 0)]
	sex <- ifelse(sex.indel == "x", "male", "female")
	change <- paste(c("is.INDEL.", "ALT.", "QUAL.", "is.INDEL.", "nb.alleles."), sex.indel, sep = "")
	for (c in change) {temp[1, c]<- NA}
	change <- paste(c("count.ref.", "count.alt.", "tot."), sex, sep = "")
	for (c in change) {temp[1, c]<- 0}
	add.INDEL.in.one.sexe <- rbind(add.INDEL.in.one.sexe, temp)
  }

  # Complete counts to 0 for variants that are found in one sexe
  add.in.one.sexe <- NULL
  for(i in 1:dim(tmp.all.SNP.in.one.sexe)[1]){
	temp <- as.data.frame(tmp.all.in.one.sexe[i, ])
        sex.no <- c("x", "y")[which(is.na(temp[, c("is.INDEL.x", "is.INDEL.y")]))]
	sex <- ifelse(sex.no == "x", "male", "female")
	change <- paste(c("count.ref.", "count.alt.", "tot."), sex, sep = "")
	for (c in change) {temp[1, c]<- 0}
	add.in.one.sexe <- rbind(add.in.one.sexe, temp)
  }
  # From these pick variants that are SNP
  ind.keep <- which(add.in.one.sexe$is.INDEL.x == 0 | add.in.one.sexe$is.INDEL.y == 0)
  add.in.one.sexe.SNP <- add.in.one.sexe[ind.keep, ]
  # From these pick variants that are INDEL
  ind.remove <- setdiff(1:dim(add.in.one.sexe)[1], ind.keep)
  add.in.one.sexe.INDEL <- add.in.one.sexe[ind.remove, ]

  ##### Final event of SNP to test for INDEL presence
  tmp.all.SNP.test.INDEL <- rbind(tmp.all.SNP, add.SNP.in.one.sexe, add.in.one.sexe.SNP)
  tmp.all.INDEL.test.INDEL <- rbind(tmp.all.INDEL, add.in.one.sexe.INDEL)

  ##### Compute distance to the nearest INDELs to remove SNP at 1bp away from INDEL and remove INDELs
 contigs <- unique(tmp.all.SNP.test.INDEL$CHROM)
  remove.id <- c()
  for (j in seq_along(contigs)){
	tmp.chr <- filter(tmp.all.SNP.test.INDEL, CHROM == contigs[j])
	tmp.chr.indel <- filter(tmp.all.INDEL.test.INDEL, CHROM == contigs[j])
	if(dim(tmp.chr.indel)[1]>0){
		dist <- apply(tmp.chr, 1, function(x) min(abs(tmp.chr.indel$POS-as.numeric(x[["POS"]]))))
		ind <- which(dist<=1)
		if(length(ind)>0){
		        remove.id <- c(remove.id, tmp.chr$ID[ind])
		}
	}
  }

  tmp.all.INDEL.test.INDEL.filt <- rbind(tmp.all.INDEL.test.INDEL, add.INDEL.in.one.sexe)
  ##### Remove SNP at 1bp away from INDEL
  tmp.all.SNP.test.INDEL.filt <- tmp.all.SNP.test.INDEL[-which(tmp.all.SNP.test.INDEL$ID %in% remove.id), ]
  ##### Add same.ALT column if ALT allele differs between male and female
  tmp.all.SNP.test.INDEL.filt$same.ALT <- 1
  tmp.all.SNP.test.INDEL.filt$same.ALT[which((tmp.all.SNP.test.INDEL.filt$ALT.x == tmp.all.SNP.test.INDEL.filt$ALT.y) == F)] <- 0
  tmp.all.INDEL.test.INDEL.filt$same.ALT <- 1
  tmp.all.INDEL.test.INDEL.filt$same.ALT[which((tmp.all.INDEL.test.INDEL.filt$ALT.x == tmp.all.INDEL.test.INDEL.filt$ALT.y) == F)] <- 0

  ##### Add contig length
  mat.size.ctg <- matrix(size.genome.mbelari$V2, nrow = 1, dimnames = list(NULL, size.genome.mbelari$V1))
  tmp.all.SNP.test.INDEL.filt$contig_length <- mat.size.ctg[,  tmp.all.SNP.test.INDEL.filt$CHROM]
  tmp.all.INDEL.test.INDEL.filt$contig_length <- mat.size.ctg[, tmp.all.INDEL.test.INDEL.filt$CHROM]

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
