require(data.table)

size.genome.mbelari <- read.table("data/ReferenceGenomes/2017_09_13_Mbelari.sizes.genome", h = F, sep = "\t", stringsAsFactors = F)

alpha <- 0.05 # treshold of significancy for adjusted pvalues

##### Quick analysis of fisher tests
data.raw.fisher <- read.csv("results/call_var/2017_11_10_merge.sexe.common.biallelic.filt.SNP.fisher_raw_count.txt", sep = "\t", h = T, stringsAsFactors = F)

# Table of significant test at ampha treshold
summary.raw.fisher <- table(data.raw.fisher$pval.adjBH.fisher.raw.count <= alpha, data.raw.fisher$CHROM)
signif.BH.raw.fisher <- which(summary.raw.fisher[2, ] >0)
nb.signif.BH.raw.fisher <- length(which(data.raw.fisher$pval.adjBH.fisher.raw.count <= alpha)) 
ctgs.signif.BH.raw.fisher <- names(summary.raw.fisher[2, signif.BH.raw.fisher])

# Nb significant tests per contig
pdf("results/call_var/2017_11_10_output_fisher_raw.SNP.fisher_raw_count.pdf", w = 10)
barplot(summary.raw.fisher[2, signif.BH.raw.fisher], cex.axis =0.4, las = 2, cex.names = 0.4, xlab = "At least one SNP signif. on this contig (BH Fisher on raw data)", ylab = "nb signif per contig", main = paste("Output of signif. fisher tests (", nb.signif.BH.raw.fisher, ") per SNP\nat a threshold of ", 100*alpha, "% per contig (", length(signif.BH.raw.fisher), ")", sep = ""))
dev.off()
size.ctgs.signif.BH.raw.fisher <-sapply(ctgs.signif.BH.raw.fisher, function(x) size.genome.mbelari[which(size.genome.mbelari$V1  == x), ]$V2) # all>1000
ctg.signif.BH.raw.fisher.max.pos <- names( which.max(summary.raw.fisher[2, signif.BH.raw.fisher]))

# pvalues along the contig with the more significant markers
output_fisher_max_pos <- data.raw.fisher[which(data.raw.fisher$CHROM == ctg.signif.BH.raw.fisher.max.pos), ]

pdf("results/call_var/2017_11_10_output_fisher_raw.SNP_highest_signif_markers_contig.pdf", w = 10)
plot(output_fisher_max_pos$POS, -log10(output_fisher_max_pos$pval.adjBH.fisher.raw.count), xlab = "Position", ylab = "-log10(pval BH fisher on raw count", main = paste("-log10(pval BH) on ",  ctg.signif.BH.raw.fisher.max.pos, "\n(marker with the highest nb of significant SNPss)", sep = ""), type = "h")
abline(h = -log10(alpha), col = "red")
dev.off()
