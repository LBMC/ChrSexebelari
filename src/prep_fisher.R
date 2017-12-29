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
  SNP.common.pos <- length(which(merge.sexe.common$is.INDEL.x == merge.sexe.common$is.INDEL.y & merge.sexe.common$is.INDEL.x == 0))
  INDEL.common.pos <- length(which(merge.sexe.common$is.INDEL.x == merge.sexe.common$is.INDEL.y & merge.sexe.common$is.INDEL.x == 1))
  SNP.common.pos.gene <- length(which(merge.sexe.common$is.INDEL.x == merge.sexe.common$is.INDEL.y & merge.sexe.common$is.INDEL.x == 0 & merge.sexe.common$genes != ""))
  INDEL.common.pos.gene <- length(which(merge.sexe.common$is.INDEL.x == merge.sexe.common$is.INDEL.y & merge.sexe.common$is.INDEL.x == 1 & merge.sexe.common$genes != ""))

  ##### All variants in both pools 
  merge.sexe.all <- merge(tab.variants.male[, c("CHROM", "POS", "REF", "ALT", "QUAL", "nb.alleles", "count.ref.male", "count.alt.male", "tot.male", "is.INDEL","genes")], tab.variants.female[, c("CHROM", "POS", "REF", "ALT", "QUAL", "nb.alleles", "count.ref.female", "count.alt.female", "tot.female", "is.INDEL", "genes")], by = c("CHROM",  "POS", "REF", "genes"), all.x = T, all.y = T)

  ##### Variants with same position and same alt allele in both pools
  merge.sexe.common.biallelic <- merge(tab.variants.male[, c("CHROM", "POS", "REF", "ALT", "QUAL", "nb.alleles", "count.ref.male", "count.alt.male", "tot.male","is.INDEL", "genes")], tab.variants.female[, c("CHROM", "POS", "REF", "ALT", "QUAL", "nb.alleles", "count.ref.female", "count.alt.female", "tot.female", "is.INDEL", "genes")], by = c("CHROM",  "POS", "REF", "ALT", "genes"))
  SNP.common.biallelic.pos <- length(which(merge.sexe.common.biallelic$is.INDEL.x == merge.sexe.common.biallelic$is.INDEL.y & merge.sexe.common.biallelic$is.INDEL.x == 0))
  INDEL.common.biallelic.pos <- length(which(merge.sexe.common.biallelic$is.INDEL.x == merge.sexe.common.biallelic$is.INDEL.y & merge.sexe.common.biallelic$is.INDEL.x == 1))
  SNP.common.biallelic.pos.gene <- length(which(merge.sexe.common.biallelic$is.INDEL.x == merge.sexe.common.biallelic$is.INDEL.y & merge.sexe.common.biallelic$is.INDEL.x == 0 & merge.sexe.common.biallelic$genes != ""))
  INDEL.common.biallelic.pos.gene <- length(which(merge.sexe.common.biallelic$is.INDEL.x == merge.sexe.common.biallelic$is.INDEL.y & merge.sexe.common.biallelic$is.INDEL.x == 1 & merge.sexe.common.biallelic$genes != ""))

  ##### Apply filter on DP (>1 in total count) and distance to INDEL:
  # for variants with same position and same alt allele in both pools
  tmp.common.biallelic <- merge.sexe.common.biallelic[which(merge.sexe.common.biallelic$tot.male> 1 | merge.sexe.common.biallelic$tot.female>1), ]
  tmp.common.biallelic$ID <- 1:dim(tmp.common.biallelic)[1]
  SNP.common.biallelic.pos.filt <- length(which(tmp.common.biallelic$is.INDEL.x == tmp.common.biallelic$is.INDEL.y & tmp.common.biallelic$is.INDEL.x == 0))
  INDEL.common.biallelic.pos.filt <- length(which(tmp.common.biallelic$is.INDEL.x == tmp.common.biallelic$is.INDEL.y & tmp.common.biallelic$is.INDEL.x == 1))
  SNP.common.biallelic.pos.gene.filt <- length(which(tmp.common.biallelic$is.INDEL.x == tmp.common.biallelic$is.INDEL.y & tmp.common.biallelic$is.INDEL.x == 0 & tmp.common.biallelic$genes != ""))
  INDEL.common.biallelic.pos.gene.filt <- length(which(tmp.common.biallelic$is.INDEL.x == tmp.common.biallelic$is.INDEL.y & tmp.common.biallelic$is.INDEL.x == 1 & tmp.common.biallelic$genes != ""))
  # for all variants in both pools
  tmp.all <- merge.sexe.all[which(merge.sexe.all$tot.male>1 | merge.sexe.all$tot.female>1), ]
  tmp.all$ID <- 1:dim(tmp.all)[1]
  SNP.common.pos.filt <- length(which(tmp.all$is.INDEL.x == tmp.all$is.INDEL.y & tmp.all$is.INDEL.x == 0))
  INDEL.common.pos.filt <- length(which(tmp.all$is.INDEL.x == tmp.all$is.INDEL.y & tmp.all$is.INDEL.x == 1))
  SNP.common.pos.gene.filt <- length(which(tmp.all$is.INDEL.x == tmp.all$is.INDEL.y & tmp.all$is.INDEL.x == 0 & tmp.all$genes != ""))
  INDEL.common.pos.gene.filt <- length(which(tmp.all$is.INDEL.x == tmp.all$is.INDEL.y & tmp.all$is.INDEL.x == 1 & tmp.all$genes != ""))

  ###### Pick variants when variants is found in one sexe (create same_class_variant = NA) or found in both but with different type like SNP in male and INDEL in female (same_class_variant = FALSE) and add this information in same_class_variant new column
  test <- apply(tmp.all[, paste("is.INDEL.", c("x", "y"), sep = "")], 1, function(x) ifelse(is.na(x[["is.INDEL.x"]]) | is.na(x[["is.INDEL.y"]]), NA, 
x[["is.INDEL.x"]] == x[["is.INDEL.y"]]))
  tab.test <- table(test, useNA = "always") 
  tmp.all$same_class_variant <- test
  tmp.all.SNP <- tmp.all[which(tmp.all$is.INDEL.x == 0 & tmp.all$is.INDEL.y == 0), ] # All filtered SNPs
  tmp.all.INDEL <- tmp.all[which(tmp.all$is.INDEL.x == 1 & tmp.all$is.INDEL.y == 1), ] # All filtered INDELs
  tmp.all.SNP.in.one.sexe <- tmp.all[which(tmp.all$same_class_variant == F), ] 
  tmp.all.in.one.sexe <- tmp.all[which(is.na(tmp.all$same_class_variant)), ] 
  ## Add supplementary SNP for variants that are both SNP and INDEL
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
  ## Add supplementary INDEL for variants that are both SNP and INDEL 
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
  ## Complete counts to 0 for variants that are found in one sexe
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

  SNP.common.pos.filt.comb.posINDEL <- dim(tmp.all.SNP.test.INDEL.filt)[1]
  INDEL.common.pos.filt.comb.posINDEL <- dim(tmp.all.INDEL.test.INDEL.filt)[1]
  SNP.common.pos.gene.filt.comb.posINDEL <- length(which(tmp.all.SNP.test.INDEL.filt$genes != ""))
  INDEL.common.pos.gene.filt.comb.posINDEL <- length(which(tmp.all.INDEL.test.INDEL.filt$genes != ""))

  ##### write summary
  summary.filt.variants <- data.frame(label = c("nb.common.pos", "nb.common.pos.gene", "nb.common.pos.same.ALT", "nb.common.pos.same.ALT.gene", "nb.common.pos.same.ALT.filt.DP>1", "nb.common.pos.same.ALT.gene.filt.DP>1", "nb.common.pos.filt.DP>1", "nb.common.pos.gene.filt.DP>1", "nb.common.pos.filt.DP>1.final", "nb.common.pos.gene.filt.DP>1.final"),
nb_SNP_INDEL = c(paste(SNP.common.pos, INDEL.common.pos, sep = "/"), paste(SNP.common.pos.gene, INDEL.common.pos.gene, sep = "/"), paste(SNP.common.biallelic.pos, INDEL.common.biallelic.pos, sep = "/"), paste(SNP.common.biallelic.pos.gene, INDEL.common.biallelic.pos.gene, sep = "/"), paste(SNP.common.biallelic.pos.filt, INDEL.common.biallelic.pos.filt, sep = "/"), paste(SNP.common.biallelic.pos.gene.filt, INDEL.common.biallelic.pos.gene.filt, sep = "/"), paste(SNP.common.pos.filt, INDEL.common.pos.filt, sep = "/"), paste(SNP.common.pos.gene.filt, INDEL.common.pos.gene.filt, sep = "/"), paste(SNP.common.pos.filt.comb.posINDEL, INDEL.common.pos.filt.comb.posINDEL, sep = "/"), paste(SNP.common.pos.gene.filt.comb.posINDEL, INDEL.common.pos.gene.filt.comb.posINDEL, sep = "/")))
pdf("results/call_var/summary_filt_variants.pdf", height = 4, width = 5)
  grid.table(summary.filt.variants)
  dev.off()
  system("bash src/date.sh results/call_var/summary_filt_variants.pdf")

  ##### Summary of filtered variants
  tab.snp <- table(tmp.all.SNP.test.INDEL.filt$same.ALT)
  tab.indel <- table(tmp.all.INDEL.test.INDEL.filt$same.ALT)
  tab.snp.g <- table(tmp.all.SNP.test.INDEL.filt[which(tmp.all.SNP.test.INDEL.filt$gene != ""), ]$same.ALT)
  tab.indel.g <- table(tmp.all.INDEL.test.INDEL.filt[which(tmp.all.INDEL.test.INDEL.filt$gene != ""), ]$same.ALT)
  pdf("results/call_var/summary_alleles_ALT_between_pools.pdf", height=8, width=10)
  par(mfrow = c(2, 2))
  bb.snp <- barplot(tab.snp, names.arg = c("different ALT allele", "same ALT allele"), main = "SNP identity between both pools\nfor all SNP", cex.main = 0.75)
  text(bb.snp, max(tab.snp)/2,labels=tab.snp) 
  bb.indel <- barplot(tab.indel, names.arg = c("different ALT allele", "same ALT allele"), main = "INDELs identity between both pools\nfor all INDEL", cex.main = 0.75) 
  text(bb.indel, max(tab.indel)/2,labels=tab.indel) 

  bb.snp.g <- barplot(tab.snp.g, names.arg = c("different ALT allele", "same ALT allele"), main = "SNP identity between both pools\nwithin genes", cex.main = 0.75)
  text(bb.snp.g, max(tab.snp.g)/2,labels=tab.snp.g) 
  bb.indel.g <- barplot(tab.indel.g, names.arg = c("different ALT allele", "same ALT allele"), main = "INDELs identity between both pools\nwithin genes", cex.main = 0.75) 
  text(bb.indel.g, max(tab.indel.g)/2,labels=tab.indel.g) 
  dev.off()
  system("bash src/date.sh results/call_var/summary_alleles_ALT_between_pools.pdf")

  write.table(tmp.all.INDEL.test.INDEL.filt, "results/call_var/merge.sexe.all.filt.INDEL.txt", sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(tmp.all.SNP.test.INDEL.filt, "results/call_var/merge.sexe.all.filt.SNP.txt", sep = "\t", row.names = F, col.names = T, quote = F)

  system("bash src/date.sh results/call_var/merge.sexe.all.filt.INDEL.txt")
  system("bash src/date.sh results/call_var/merge.sexe.all.filt.SNP.txt")
}

##### Summary of density of SNPs for filtered SNPs per contig 
tests  <- fread("results/call_var/2017_12_08_merge.sexe.all.filt.SNP.txt", sep = "\t", h = T)

# At contig
tmp.tests <- tests[which(tests$tot.female>0), ]
nb.test.female <- table(tmp.tests$CHROM)
tmp.tests <- tests[which(tests$tot.male>0), ]
nb.test.male <- table(tmp.tests$CHROM)
density.tests.sexe.ctg <- data.frame(contig = size.genome.mbelari$V1, nb.female= as.vector(unname(nb.test.female[size.genome.mbelari$V1])), nb.male = as.vector(unname(nb.test.male[size.genome.mbelari$V1])), stringsAsFactors = F)
density.tests.sexe.ctg$log2.den.female.male <- log2(density.tests.sexe.ctg$nb.female/density.tests.sexe.ctg$nb.male)
write.table(density.tests.sexe.ctg, "results/call_var/density.contigs.merge.sexe.all.filt.SNP.txt", sep = "\t", row.names = F, col.names = T, quote = F)
system("bash src/date.sh results/call_var/density.contigs.merge.sexe.all.filt.SNP.txt")

# At gene
g <- unique(tests$genes)
tmp.tests <- tests[which(tests$tot.female>0 & tests$genes != ""), ]
nb.test.female <- table(tmp.tests$genes)
tmp.tests <- tests[which(tests$tot.male>0  & tests$genes != ""), ]
nb.test.male <- table(tmp.tests$genes)
density.tests.sexe.g <- data.frame(gene = g, nb.female= as.vector(unname(nb.test.female[g])), nb.male = as.vector(unname(nb.test.male[g])), stringsAsFactors = F)
density.tests.sexe.g$log2.den.female.male <- log2(density.tests.sexe.g$nb.female/density.tests.sexe.g$nb.male)
write.table(density.tests.sexe.g, "results/call_var/density.genes.merge.sexe.all.filt.SNP.txt", sep = "\t", row.names = F, col.names = T, quote = F)
system("bash src/date.sh results/call_var/density.genes.merge.sexe.all.filt.SNP.txt")

# tab at ctg
pos <- density.tests.sexe.ctg$log2.den.female.male[which(density.tests.sexe.ctg$log2.den.female.male>0 & is.na(density.tests.sexe.ctg$log2.den.female.male) == F)]
neg <- density.tests.sexe.ctg$log2.den.female.male[which(density.tests.sexe.ctg$log2.den.female.male<0 & is.na(density.tests.sexe.ctg$log2.den.female.male) == F)]
data.pos <- transform(pos, groupdata = cut(pos, breaks=c(0, 0.5, 1, 1.5, 2, 2.5),
right=F, include.lowest = F))
data.neg <- transform(neg, groupdata = cut(neg, breaks=rev(-c(0, 0.5, 1, 1.5, 2, 2.5)),
right=T, include.lowest = F))
data <- rbind(data.neg, data.pos)
tab <- table(c(as.character(data$groupdata), rep("NA", length(which(is.na(density.tests.sexe.ctg$log2.den.female.male)))), rep("0", length(which(density.tests.sexe.ctg$log2.den.female.male == 0)))))
names <- c("NA", "(-2.5,-2]", "(-2,-1.5]", "(-1.5,-1]", "(-1,-0.5]", "(-0.5,0]", "0", "[0,0.5)", "[0.5,1)", "[1,1.5)", "[1.5,2)" , "[2,2.5)")
tab.ctg <- tab[names]
ind.na <- which(is.na(tab.ctg))
tab.ctg[ind.na] <- 0 
names(tab.ctg)[ind.na] <- paste(names[ind.na], "\n0 count")
names(tab.ctg)[which(names(tab.ctg) == "NA")] <- paste("0 in male\n(", tab.ctg["NA"], ")", sep = "")
names(tab.ctg)[which(names(tab.ctg) == "0")] <- "male\n=fem"

# tab at gene
pos <- density.tests.sexe.g$log2.den.female.male[which(density.tests.sexe.g$log2.den.female.male>0 & is.na(density.tests.sexe.g$log2.den.female.male) == F)]
neg <- density.tests.sexe.g$log2.den.female.male[which(density.tests.sexe.g$log2.den.female.male<0 & is.na(density.tests.sexe.g$log2.den.female.male) == F)]
data.pos <- transform(pos, groupdata = cut(pos, breaks=c(0, 0.5, 1, 1.5, 2, 2.5),
right=F, include.lowest = F))
data.neg <- transform(neg, groupdata = cut(neg, breaks=rev(-c(0, 0.5, 1, 1.5, 2, 2.5)),
right=T, include.lowest = F)) 
data <- rbind(data.neg, data.pos)
tab <- table(c(as.character(data$groupdata), rep("NA", length(which(is.na(density.tests.sexe.g$log2.den.female.male)))), rep("0", length(which(density.tests.sexe.g$log2.den.female.male == 0)))))
names <- c("NA", "(-2.5,-2]", "(-2,-1.5]", "(-1.5,-1]", "(-1,-0.5]", "(-0.5,0]", "0", "[0,0.5)", "[0.5,1)", "[1,1.5)", "[1.5,2)" , "[2,2.5)")
tab.g <- tab[names]
ind.na <- which(is.na(tab.g))
tab.g[ind.na] <- 0 
names(tab.g)[ind.na] <- paste(names[ind.na], "\n0 count")
names(tab.g)[which(names(tab.g) == "NA")] <- paste("0 in male\n(", tab.g["NA"], ")", sep = "")
names(tab.g)[which(names(tab.g) == "0")] <- "male\n=fem"

##### Summary of density of SNPs for filtered variants per contig
pdf("results/call_var/log2_density_ratio_filt_SNPs.pdf", w = 13, h = 5)
par(mfrow = c(1,2))
barplot(tab.ctg, main = "log2(density SNP female/density SNP male)\non filtered detected SNPs on contig", xlab  = "log2(density SNP female/density SNP male)", cex.names = 0.4)
barplot(tab.g, main = "log2(density SNP female/density SNP male)\non filtered detected SNPs on gene", xlab  = "log2(density SNP female/density SNP male)", cex.names = 0.4)
dev.off()
system("bash src/date.sh results/call_var/log2_density_ratio_SNPs.pdf")

##### Plot histogram of frequency of filtered variants 
pdf("results/call_var/frequency_SNPs_sexes.pdf")
mean.m <- mean(na.omit(tests$count.alt.male/tests$tot.male))
mean.f <-mean(na.omit(tests$count.alt.female/tests$tot.female))
hist(tests$count.alt.male/tests$tot.male, col = adjustcolor("blue", 0.5), breaks = 20, xlim = c(0,1), xlab = "Variants frequency", main = "Histogram of SNPs frequencies\nbetween sexes")
hist(tests$count.alt.female/tests$tot.female, col = adjustcolor("red", 0.5), add = T, breaks = 20, xlim = c(0,1))
abline(v = mean.f, col = "red", lty = 2, lwd = 2)
abline(v = mean.m, col = "blue", lty = 2, lwd = 2)
legend("topleft", c(paste("male frequency (SNPs found in\nfemale only = ", length(which(tests$tot.female == 0)), ")", sep = ""), paste("female frequency (SNPs found in\nmale only = ", length(which(tests$tot.male == 0)), ")", sep = "")), fill = c("blue", "red"))
dev.off()
system("bash src/date.sh results/call_var/frequency_SNPs_sexes.pdf")
