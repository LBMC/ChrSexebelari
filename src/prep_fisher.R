require(data.table)
require(dplyr)
library(foreach)
library(doParallel)
library(gridExtra)

size.genome.mbelari <- read.table("data/ReferenceGenomes/2017_09_13_Mbelari.sizes.genome", h = F, sep = "\t", stringsAsFactors = F)

do.prep.fisher <- F

tab.variants.female <- fread("results/call_var/2017_11_09_female_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts_genic_information.txt", sep = "\t", h = T, stringsAsFactors = F)
tab.variants.male <- fread("results/call_var/2017_11_09_male_trim_Mbelari_mapped_rmdup_rg_realign_indels_counts_genic_information.txt", sep = "\t", h = T, stringsAsFactors = F)

if (do.prep.fisher) {
 
  ##### Look at common / biallelic genic sites between pools
  merge.sexe.common <- merge(tab.variants.male[, c("CHROM", "POS", "REF", "ALT", "QUAL", "nb.alleles", "count.ref.male", "count.alt.male", "tot.male", "genes")], tab.variants.female[, c("CHROM", "POS", "REF", "ALT", "QUAL", "nb.alleles", "count.ref.female", "count.alt.female", "tot.female", "is.INDEL", "genes")], by = c("CHROM",  "POS", "REF", "genes"))
  ind.genic <- which(merge.sexe.common$genes != "")
  nb.SNPs.all.common = table(merge.sexe.common$is.INDEL)["0"]
  nb.INDELs.all.common =  table(merge.sexe.common$is.INDEL)["1"]
  nb.SNPs.genic.common = table(merge.sexe.common[ind.genic]$is.INDEL)["0"]
  nb.INDELs.genic.common = table(merge.sexe.common[ind.genic]$is.INDEL)["1"]
  nb.ctg.SNPs.all.common <- length(unique(merge.sexe.common$CHROM[which(merge.sexe.common$is.INDEL == 0)]))
  nb.ctg.INDELs.all.common <- length(unique(merge.sexe.common$CHROM[which(merge.sexe.common$is.INDEL == 1)]))

  merge.sexe.common.biallelic <- merge(tab.variants.male[, c("CHROM", "POS", "REF", "ALT", "QUAL", "nb.alleles", "count.ref.male", "count.alt.male", "tot.male", "genes")], tab.variants.female[, c("CHROM", "POS", "REF", "ALT", "QUAL", "nb.alleles", "count.ref.female", "count.alt.female", "tot.female", "is.INDEL", "genes")], by = c("CHROM",  "POS", "REF", "ALT", "genes"))
  ind.genic <- which(merge.sexe.common.biallelic$genes != "")
  nb.SNPs.all.common.biallelic = table(merge.sexe.common.biallelic$is.INDEL)["0"]
  nb.INDELs.all.common.biallelic =  table(merge.sexe.common.biallelic$is.INDEL)["1"]
  nb.SNPs.genic.common.biallelic = table(merge.sexe.common.biallelic[ind.genic]$is.INDEL)["0"]
  nb.INDELs.genic.common.biallelic = table(merge.sexe.common.biallelic[ind.genic]$is.INDEL)["1"]
  nb.ctg.SNPs.all.common.biallelic <- length(unique(merge.sexe.common.biallelic$CHROM[which(merge.sexe.common.biallelic$is.INDEL == 0)]))
  nb.ctg.INDELs.all.common.biallelic <- length(unique(merge.sexe.common.biallelic$CHROM[which(merge.sexe.common.biallelic$is.INDEL == 1)]))

  ##### Apply filter on contig size, DP (>1 in each pool) and distance to INDEL
  tmp <- merge.sexe.common.biallelic[which(merge.sexe.common.biallelic$tot.male>1 & merge.sexe.common.biallelic$tot.female>1), ]
  tmp$ID <- 1:dim(tmp)[1]
	
  size.contigs1000 <- size.genome.mbelari$V1[which(size.genome.mbelari$V2 >= 1000)]
  tmp <- tmp[which(tmp$CHROM %in% size.contigs1000), ]

  ##### Compute distance to the nearest INDELs to remove SNP at 1bp away from INDEL and remove INDELs
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  remove.dist <- foreach(i=1:dim(tmp)[1], .combine=c) %dopar% {
  	tmp2 <- tmp[i, ]
  	if(tmp2$is.INDEL == 0) {
		temp <- tmp[which(tmp$CHROM == tmp2$CHROM & tmp$is.INDEL == 1), ]
		d <- temp$POS-tmp2$POS
		dist <- min(abs(d))
		if (dist <= 1){
			remove <- tmp2$ID
		}else{
			remove <- NA
		}
	} else{
		remove <- NA
	}
	remove
  }
  stopCluster(cl)

  ids.remove <- remove.dist[-which(is.na(remove.dist))]
  merge.sexe.common.biallelic.filt <- tmp[-which(tmp$ID %in% ids.remove), ]
  ind.genic <- which(merge.sexe.common.biallelic.filt$genes != "")
  nb.SNPs.all.common.biallelic.filt = table(merge.sexe.common.biallelic.filt$is.INDEL)["0"]
  nb.INDELs.all.common.biallelic.filt =  table(merge.sexe.common.biallelic.filt$is.INDEL)["1"]
  nb.SNPs.genic.common.biallelic.filt = table(merge.sexe.common.biallelic.filt[ind.genic]$is.INDEL)["0"]
  nb.INDELs.genic.common.biallelic.filt = table(merge.sexe.common.biallelic.filt[ind.genic]$is.INDEL)["1"]
  nb.ctg.SNPs.all.common.biallelic.filt <- length(unique(merge.sexe.common.biallelic.filt$CHROM[which(merge.sexe.common.biallelic.filt$is.INDEL == 0)]))
  nb.ctg.INDELs.all.common.biallelic.filt <- length(unique(merge.sexe.common.biallelic.filt$CHROM[which(merge.sexe.common.biallelic.filt$is.INDEL == 1)]))

  summary.filt <- data.frame(nb.SNPs.all.common = nb.SNPs.all.common, nb.INDELs.all.common =  nb.INDELs.all.common, nb.SNPs.genic.common = nb.SNPs.genic.common,
  nb.INDELs.genic.common = nb.INDELs.genic.common, nb.ctg.SNPs.all.common = nb.ctg.SNPs.all.common, nb.ctg.INDELs.all.common = nb.ctg.INDELs.all.common,
  nb.SNPs.all.common.biallelic = nb.SNPs.all.common.biallelic, nb.INDELs.all.common.biallelic =  nb.INDELs.all.common.biallelic, nb.SNPs.genic.common.biallelic = nb.SNPs.genic.common.biallelic, nb.INDELs.genic.common.biallelic = nb.INDELs.genic.common.biallelic, nb.ctg.SNPs.all.common.biallelic = nb.ctg.SNPs.all.common.biallelic,
  nb.ctg.INDELs.all.common.biallelic = nb.ctg.INDELs.all.common.biallelic, nb.SNPs.all.common.biallelic.filt = nb.SNPs.all.common.biallelic.filt,
  nb.INDELs.all.common.biallelic.filt =  nb.INDELs.all.common.biallelic.filt, nb.SNPs.genic.common.biallelic.filt = nb.SNPs.genic.common.biallelic.filt,
  nb.INDELs.genic.common.biallelic.filt = nb.INDELs.genic.common.biallelic.filt, nb.ctg.SNPs.all.common.biallelic.filt = nb.ctg.SNPs.all.common.biallelic.filt,
  nb.ctg.INDELs.all.common.biallelic.filt = nb.ctg.INDELs.all.common.biallelic.filt, stringsAsFactors = F)

  write.table(summary.filt, "results/call_var/2017_11_10_summary_filt_variants.txt", sep = "\t",  quote = F, row.names = F)
  write.table(merge.sexe.common.biallelic.filt, "results/call_var/2017_11_09_merge.sexe.common.biallelic.filt.txt", sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(merge.sexe.common.biallelic.filt[which(merge.sexe.common.biallelic.filt$is.INDEL == 0),], "results/call_var/2017_11_09_merge.sexe.common.biallelic.filt.SNP.txt", sep = "\t", row.names = F, col.names = T, quote = F)
}

f <- read.table("results/call_var/2017_11_10_summary_filt_variants.txt", sep = "\t", h = T)
pdf("results/call_var/2017_11_10_summary_filt1_variants.pdf", height=4, width=14)
grid.table(f[, 1:6])
dev.off()

pdf("results/call_var/2017_11_10_summary_filt2_variants.pdf", height=4, width=14)
grid.table(f[, 7:12])
dev.off()

pdf("results/call_var/2017_11_10_summary_filt3_variants.pdf", height=4, width=14)
grid.table(f[, 13:18])
dev.off()
