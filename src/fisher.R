library(gridExtra)
require(modeest)
require(Sushi)
require(qvalue)
require(data.table)
require(foreach)
require(doParallel)
require(VGAM)
source("src/func/functions.R")

do.fisher <- F
alpha <- 0.05

if (do.fisher) {
  data.filt.SNP <- read.csv("results/call_var/2017_12_08_merge.sexe.all.filt.SNP.txt", sep = "\t", h = T, stringsAsFactors = F)

  # Fisher tests are not computed on all SNPs!!
  # Pick same alt allele for female and male and with contig length >= 1000
  data.filt.SNP.same.alt <- data.filt.SNP[which(data.filt.SNP$same.ALT == 1 & data.filt.SNP$contig_length>=1000), ]

  # Compute freq male and freq female
  data.filt.SNP.same.alt$freq.alt.female <- data.filt.SNP.same.alt$count.alt.female/(data.filt.SNP.same.alt$count.ref.female+data.filt.SNP.same.alt$count.alt.female)
  data.filt.SNP.same.alt$freq.alt.male <- data.filt.SNP.same.alt$count.alt.male/(data.filt.SNP.same.alt$count.ref.male+data.filt.SNP.same.alt$count.alt.male)

  ##################
  ##### Fisher #####
  ##################

  # Add raw pvalue associated to Fisher test done on raw counts
  pval.fisher <- ComputePvalueFisher(input.df = data.filt.SNP.same.alt[, c("count.ref.male", "count.alt.male", "count.ref.female", "count.alt.female")])
  data.filt.SNP.same.alt$pval.fisher.raw.count <- pval.fisher

  # Adjust Fisher pvalue with BH
  data.filt.SNP.same.alt$pval.adjBH.fisher.raw.count <- p.adjust(data.filt.SNP.same.alt$pval.fisher.raw.count, method = "BH")

  write.table(data.filt.SNP.same.alt, "results/call_var/merge.sexe.all.filt.SNP.same.alt.fisher.txt", sep = "\t",  quote = F, row.names = F)
  system("bash src/date.sh results/call_var/merge.sexe.all.filt.SNP.same.alt.fisher.txt")
}

do.H0 <- F
if(do.H0){
	##### Compute H0 distribution of count tables: UPDATE: instead using Bastide et al Methods to choose the best alpha value, I use instead a vglm fit on betabinomial to get alpha and beta parameters.
	data.tests <- read.csv("results/call_var/2017_12_05_merge.sexe.common.biallelic.filt.SNP.same.alt.fisher.txt", sep = "\t", h = T, stringsAsFactors = F)
	fit <- vglm(cbind(c(data.tests$count.ref.female, data.tests$count.alt.female), c(data.tests$count.ref.male, data.tests$count.alt.male)) ~ 1, betabinomial, trace = T)
	coeff <- Coef(fit)
	alpha <- (1-coeff["rho"])/(coeff["rho"]*(1+(1-coeff["mu"])/coeff["mu"]))
	beta <- alpha*(1-coeff["mu"])/coeff["mu"]
	saveRDS(fit, "results/call_var/params_H0_beta.RData")
	system("bash src/date.sh results/call_var/params_H0_beta.RData")
}

##### Compute 10 datasets under H0 with alpha set to obtain above
do.simu.10x <- F
if(do.simu.10x) {
	data.tests <- read.csv("results/call_var/2017_12_05_merge.sexe.common.all.filt.SNP.same.alt.fisher.txt", sep = "\t", h = T, stringsAsFactors = F)
        ##### complete
	alpha.choose <- alpha
	beta.choose <- beta
	ind.ord <- order(data.tests$pval.fisher.raw.count)
	data.tests.simu <- data.tests[ind.ord, ]
	pH0.simu <- NULL
	for (i in 1:10){
		print(i)
		cores=detectCores()
		cl <- makeCluster(cores[1]-1) #not to overload your computer
		registerDoParallel(cl)
		pvalH0_out <- foreach(i=1:dim(data.tests.simu)[1], .combine=rbind) %dopar% {
	          	tmp <- data.tests.simu[i, ]
	
			nref <- tmp$count.ref.female + tmp$count.ref.male
	                pmref <- rbeta(1, shape1 = alpha.choose, shape2 = beta.choose)
			nrefmale <- rbinom(n = 1, size = nref, prob = pmref)
			nreffemale <- nref-nrefmale
	
			nalt <- tmp$count.alt.female + tmp$count.alt.male
                	pmalt <- rbeta(1, shape1 = alpha.choose, shape2 = beta.choose)
			naltmale <- rbinom(n = 1, size = nalt, prob = pmalt)
			naltfemale <- nalt-naltmale

 			data.frame(ID = tmp$ID, count.ref.male = nrefmale, count.alt.male = naltmale, count.ref.female = nreffemale, count.alt.female = naltfemale, stringsAsFactors = F)
		}
	  	stopCluster(cl)

		pvalH0_out <- pvalH0_out[sapply(data.tests.simu$ID, function(x) which(pvalH0_out$ID == x)), ]
	        pval.fisher.H0 <- ComputePvalueFisher(input.df = pvalH0_out[, c("count.ref.male", "count.alt.male", "count.ref.female", "count.alt.female")])
	        pvalH0_out$pval.fisher.H0 <- pval.fisher.H0
		pH0.simu <- cbind(pH0.simu, pval.fisher.H0)

		saveRDS(pvalH0_out, paste("results/fisher/simulated_#", i, "_count_H0.RData", sep = ""))	
		system(paste("bash src/date.sh results/fisher/simulated_#", i, "_count_H0.RData", sep = ""))	
	}

	pH0.simu <- cbind(pH0.simu, data.frame(pval.obs.ord = data.tests.simu$pval.fisher.raw.count))
	# For the ith rank SNP, # nb simu SNP with pvalue <= pvalue obs for Fisher and compute fdr with rank
	pemp <- apply(pH0.simu, 1, function(x) {y <- unlist(x); length(which(y[1:10] <= y[11]))})/10
	fdr <- p.adjust(pemp, method = "BH")
	data.tests.simu$FDR <- fdr
	write.table(data.tests.simu, "results/call_var/merge.sexe.common.all.filt.SNP.same.alt.fisher.FDR.txt", sep = "\t", col.names = T, row.names = F)
	system("bash src/date.sh results/call_var/merge.sexe.common.all.filt.SNP.same.alt.fisher.FDR.txt")
}	

##### Zoom on tests output
size.genome.mbelari <- fread("data/ReferenceGenomes/2017_09_13_Mbelari.sizes.genome", sep = "\t", h = F, stringsAsFactors = F)
# annotation
gff.mbelari <- read.csv("data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2.gff3", h = F, sep = "\t", stringsAsFactors = F)
# genes annotation
gff.mbelari.genes <- gff.mbelari[which(gff.mbelari$V3 == "gene"), ]
tests <- read.csv("results/call_var/2017_12_05_merge.sexe.common.all.filt.SNP.same.alt.fisher.FDR.txt", sep = "\t", h = T, stringsAsFactors = F)

##### Nb of significative SNP results at gene and contig level
tab.ctg <- table(tests$FDR<alpha, tests$CHROM); tab.signif.ctg <- tab.ctg[, which(tab.ctg[2, ]>0)]
tab.gene <- table(tests$FDR<alpha, tests$gene); tab.signif.gene <- tab.gene[, which(tab.gene[2, ]>0)]
pdf("results/call_var/signif_Fisher_contig_gene_level.pdf")
par(mfrow = c(1,2))
barplot(tab.signif.ctg[2, ], xlab = "At least 1 SNP signif. on this contig (FDR for fisher)", ylab = "#signif per contig", xaxt = "n",
main = paste("Output of signif. fisher tests (", sum(tab.signif.ctg[2, ]), ") per\nSNP at a threshold of ", 100*alpha, "% per contig (", dim(tab.signif.ctg)[2], ")", sep = ""))
barplot(tab.signif.gene[2, ], xlab = "At least 1 SNP signif. on this gene (FDR for fisher)", ylab = "#signif per gene", xaxt = "n",
main = paste("Output of signif. fisher tests (", sum(tab.signif.gene[2, ]), ") per\nSNP at a threshold of ", 100*alpha, "% per gene (", dim(tab.signif.gene)[2], ")", sep = ""))
dev.off()
system("bash src/date.sh results/call_var/signif_Fisher.pdf")

##### Plot FDR representation by contig with highest nb of significant SNPs to lowest number of significant SNPs
fisher <- tests[, c("CHROM", "POS", "POS", "REF", "pval.fisher.raw.count", "pval.adjBH.fisher.raw.count", "FDR")]
tab <- table(tests$FDR<alpha, tests$CHROM)
tab <- tab[1:2, which(tab[2, ]>0)]
length.tab <- sapply(colnames(tab), function(x) size.genome.mbelari[which(size.genome.mbelari$V1 == x), ]$V2)
nb.per.bp <- tab[2, ]/length.tab[which(tab[2, ]>0)]
nb.on.tot <-100*tab[2, ]/(tab[2,]+tab[1,])

##### Ranking on contigs
rank.signif.ctgs <- colnames(tab[1:2, rev(order(tab[2, ]))])
rank.signif.ctgs.per.bp <- names(rev(order(nb.per.bp)))
rank.signif.ctgs.on.tot <- names(rev(order(nb.on.tot)))

nb.signif <- 10
to.plot <- rank.signif.ctgs[1:nb.signif]
to.plot.per.bp <- rank.signif.ctgs.per.bp[1:nb.signif]
to.plot.on.tot <- rank.signif.ctgs.on.tot[1:nb.signif]
plot <- to.plot

# Prepare input for Manhattan plot
size.tmp <- size.genome.mbelari[sapply(plot, function(x) which(size.genome.mbelari$V1 == x)), ]
size.tmp$V1 <- factor(size.tmp$V1, levels =  plot)
fisher.tmp <- fisher[unlist(sapply(plot, function(x) which(fisher$CHROM == x))), c("CHROM", "POS", "POS.1", "REF", "pval.adjBH.fisher.raw.count")]
fisher.tmp$CHROM <- paste("chr", fisher.tmp$CHROM, sep = "")
fisher.tmp$CHROM <- factor(fisher.tmp$CHROM, levels = paste("chr", plot, sep = ""))
fisher.tmp$REF <- 1:dim(fisher.tmp)[1]
#levels(fisher.tmp$CHROM) <- plot

pdf("results/call_var/2017_11_26_example_Manhattan.pdf")
plotManhattan(bedfile = fisher.tmp, pvalues = fisher.tmp$FDR,
col = SushiColors(5),  genome = size.tmp, cex=0.75)
labelgenome(genome = size.tmp, n = 5, edgeblankfraction = 0.20, cex.axis=.5)
abline(h = -log10(0.05))
axis(side=2,las=2,tcl=.2)
mtext("-log10(P)",side=2,line=1.75,cex=1,font=2)
mtext("contig",side=1,line=1.75,cex=1,font=2)
dev.off()
