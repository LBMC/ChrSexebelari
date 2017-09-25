require(data.table)
#require(qvalue)

size.genome.mbelari <- read.table("data/ReferenceGenomes/2017_09_13_Mbelari.sizes.genome", h = F, sep = "\t", stringsAsFactors = F)
do.computation.per.pool <- F

if (do.computation.per.pool) {
  # Unzip vcf files
  system("gunzip -c results/call_var/2017_09_21_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz > results/call_var/2017_09_21_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf")
  system("gunzip -c results/call_var/2017_09_21_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf.gz > results/call_var/2017_09_21_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf")
  
  variants.female = fread("results/call_var/2017_09_21_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf", sep = "\t",  sep2 = ";", h = T, skip = 19286)
  variants.male = fread("results/call_var/2017_09_21_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels.vcf", sep = "\t",  sep2 = ";", h = T, skip = 19286)

  # Get nb of alleles per pool (AN), number of reads "equal" to ref/number of reads equal to alt at a given position (DP4 field)
  info.female <- ExtractCountFromVCF(variants.female, "female")
  info.male <- ExtractCountFromVCF(variants.male, "male")
  
  variants.female <- cbind(variants.female, info.female)
  variants.male  <- cbind(variants.male, info.male)
  colnames(variants.female)[1] <- "CHROM"
  colnames(variants.male)[1] <- "CHROM"
  
  write.table(variants.female, "results/call_var/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts.txt", sep = "\t",  quote = F, row.names = F)
  write.table(variants.male, "results/call_var/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts.txt", sep = "\t",  quote = F, row.names = F)
}

tab.variants.female <- fread("results/call_var/MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts.txt", sep = "\t", h = T, stringsAsFactors = F)
tab.variants.male <- fread("results/call_var/MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts.txt", sep = "\t", h = T, stringsAsFactors = F)

# First work on biallelic common sites between pools
merge.sexe <- merge(variants.male[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "nb.alleles", "count.ref.male", "count.alt.male", "tot.male")], variants.female[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "nb.alleles", "count.ref.female", "count.alt.female", "tot.female", "is.INDEL")], by = c("CHROM",  "POS", "REF", "ALT", "nb.alleles"))

size.contigs1000 <- size.genome.mbelari$V1[which(size.genome.mbelari$V2 >= 1000)]
merge.sexe.filt <- merge.sexe[which(merge.sexe$CHROM %in% size.contigs1000), ]

#> dim(merge.sexe.filt)
#[1] 586769     16
#> dim(merge.sexe)
#[1] 590459     16

# Compute normalized counts by median number of counts and 1000 individuals at a given position (DP field) 
median.tot.male <- median(merge.sexe.filt$tot.male)
median.tot.female <- median(merge.sexe.filt$tot.female)

merge.sexe.filt$count.norm.ref.male <- (merge.sexe.filt$count.ref.male/median.tot.male)
merge.sexe.filt$count.norm.alt.male <- (merge.sexe.filt$count.alt.male/median.tot.male)
p.ref.male <- merge.sexe.filt$count.norm.ref.male/(merge.sexe.filt$count.norm.ref.male+merge.sexe.filt$count.norm.alt.male)
merge.sexe.filt$count.norm.ref.male <- round(1000*p.ref.male)
merge.sexe.filt$count.norm.alt.male <- 1000-merge.sexe.filt$count.norm.ref.male

merge.sexe.filt$count.norm.ref.female <- (merge.sexe.filt$count.ref.female/median.tot.female)
merge.sexe.filt$count.norm.alt.female <- (merge.sexe.filt$count.alt.female/median.tot.female)
p.ref.female <- merge.sexe.filt$count.norm.ref.female/(merge.sexe.filt$count.norm.ref.female+merge.sexe.filt$count.norm.alt.female)
merge.sexe.filt$count.norm.ref.female <- round(1000*p.ref.female)
merge.sexe.filt$count.norm.alt.female <- 1000-merge.sexe.filt$count.norm.ref.female

# Add raw pvalue associated to Fisher test done on raw counts
pval.fisher <- apply(merge.sexe.filt[, c("count.ref.male", "count.alt.male", "count.ref.female", "count.alt.female")], 1, function(x) {
f = fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol = 2, byrow = T));
pval = f$p.value;
})
merge.sexe.filt$pval.fisher.raw.count <- pval.fisher

pval.fisher.norm <- apply(merge.sexe.filt[, c("count.norm.ref.male", "count.norm.alt.male", "count.norm.ref.female", "count.norm.alt.female")], 1, function(x) {
f = fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol = 2, byrow = T));
pval = f$p.value;
})
merge.sexe.filt$pval.fisher.norm.count <- pval.fisher.norm

# Compute freq male and freq female
merge.sexe.filt$freq.alt.female <- merge.sexe.filt$count.alt.female/(merge.sexe.filt$count.ref.female+merge.sexe.filt$count.alt.female)
merge.sexe.filt$freq.alt.male <- merge.sexe.filt$count.alt.male/(merge.sexe.filt$count.ref.male+merge.sexe.filt$count.alt.male)

merge.sexe.filt$freq.norm.alt.female <- merge.sexe.filt$count.norm.alt.female/(merge.sexe.filt$count.norm.ref.female+merge.sexe.filt$count.norm.alt.female)
merge.sexe.filt$freq.norm.alt.male <- merge.sexe.filt$count.norm.alt.male/(merge.sexe.filt$count.norm.ref.male+merge.sexe.filt$count.norm.alt.male)

# Adjust pvalue with BH
merge.sexe.filt$pval.adjBH.fisher.raw.count <- p.adjust(merge.sexe.filt$pval.fisher.raw.count, method = "BH")
merge.sexe.filt$pval.adjBH.fisher.norm.count <- p.adjust(merge.sexe.filt$pval.fisher.norm.count, method = "BH")
# Compute qvalue 
#qobj.fisher <- qvalue(p = merge.sexe.filt$pval.fisher.raw.count)
#merge.sexe.filt$qval.fisher.raw.count <- qobj.fisher$qvalues
#qobj.norm.fisher <- qvalue(p = merge.sexe.filt$pval.fisher.norm.count)
#merge.sexe.filt$qval.fisher.norm.count <- qobj.norm.fisher$qvalues


# Compute average variants density as  the nb of variants per contigs divided by the number of bp in the contig with a DP !=0 

# Annotate masked regions 

# Summary of significative position by contig at a rough threshold of 0.05
contig.tests <- unique(merge.sexe.filt[, "CHROM"])
summary.tests <- t(sapply(seq_along(contig.tests), function(x) {c <- contig.tests[x];
s <- size.genome.mbelari[which(size.genome.mbelari$V1 == c), ]$V2;  
ind.contig <- which(merge.sexe.filt[, "CHROM"] == c);
c(c, s, length(ind.contig), length(which(merge.sexe.filt$pval.adjBH.fisher.raw.count[ind.contig] < 0.05)), length(which(merge.sexe.filt$pval.adjBH.fisher.norm.count[ind.contig] < 0.05)))
}))
summary.tests <- as.data.frame(summary.tests, stringsAsFactors = F)
colnames(summary.tests) <- c("contig", "length", "nb.tests", "nb.signif.tests.raw", "nb.signif.tests.norm")
summary.tests$nb.tests <- as.numeric(summary.tests$nb.tests)
summary.tests$nb.signif.tests.raw <- as.numeric(summary.tests$nb.signif.tests.raw)
summary.tests$nb.signif.tests.norm <- as.numeric(summary.tests$nb.signif.tests.norm)
summary.tests$length <- as.numeric(summary.tests$length)

ind.ord <- order(summary.tests$length, decreasing = T)#/summary.pval.norm.fisher.0.05.bonf.tests$length)
contig.ord <- summary.tests$contig[ind.ord]
contig.length.ord <- summary.tests$length[ind.ord]

# Reorder contigs file directly if you want to plot no in the same order as in the txt
which.contigs <- 1:1000
offset <- 0; xlim <- c(0, sum(contig.length.ord[which.contigs])); ylim = c(-30,30)
for (i in which.contigs) {
	c <- contig.ord[i]	
	s <- contig.length.ord[i]	
	ind <- which(merge.sexe.filt[, "#CHROM"] == c)
	tmp <- merge.sexe.filt[ind, ]
	col <- as.numeric(tmp$freq.alt.female>tmp$freq.alt.male)
	pos <- tmp$POS + offset
	offset <- offset + s
	if (i == 1) {
		#plot(pos, -log10(tmp$pval.fisher.norm.count)*c(1,-1)[col+1], col = c("blue", "red")[col+1], type = "h", xlim = xlim, ylim = ylim)
		plot(pos, -log10(tmp$pval.adjBH.fisher.norm.count)*c(1,-1)[col+1], col = c("blue", "red")[col+1], type = "p", pch = 20, xlim = xlim, ylim =  ylim)
		abline(h = 0, col = "black")
	}else{
		#points(pos, -log10(tmp$pval.fisher.norm.count)*c(1,-1)[col+1], col = c("blue", "red")[col+1], type = "h")
		points(pos, -log10(tmp$pval.adjBH.fisher.norm.count)*c(1,-1)[col+1], col = c("blue", "red")[col+1], type = "p", pch = 20)
	}
	abline(v = offset, lty = 2, col = "gray")
}

# Summary nb of tests per kb
nb.tests.per.kb <- summary.tests$nb.tests*1000/summary.tests$length
#> summary(nb.tests.per.kb)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00546  2.05931  3.56061  3.89762  5.18996 36.82171 

# Plot now log2(tot reads female/tot reads male)
which.contigs <- 1:10
offset <- 0; xlim <- c(0, sum(contig.length.ord[which.contigs])); ylim = c(-6,6)
for (i in which.contigs) {
	c <- contig.ord[i]	
	s <- contig.length.ord[i]	
	ind <- which(merge.sexe.filt[, "#CHROM"] == c)
	tmp <- merge.sexe.filt[ind, ]
	pos <- tmp$POS + offset
	offset <- offset + s
	if (i == 1) {
		plot(pos, log2((tmp$tot.female/median.tot.female)/(tmp$tot.male/median.tot.male)), type = "p", col = adjustcolor("black", 0.5), pch = 20, xlim = xlim, ylim =  ylim)
		abline(h = 0, col = "black")
	}else{
		points(pos, log2((tmp$tot.female/median.tot.female)/(tmp$tot.male/median.tot.male)), type = "p", col = adjustcolor("black", 0.5), pch = 20)
	}
	abline(v = offset, lty = 2, col = "gray")
}


# Do the same on corrected counts and with G stats

