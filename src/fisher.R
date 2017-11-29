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
  data.filt.SNP <- read.csv("results/call_var/2017_11_13_merge.sexe.all.filt.SNP.txt", sep = "\t", h = T, stringsAsFactors = F)

  # Pick same alt allele for female and male
  data.filt.SNP.same.alt <- data.filt.SNP[which(data.filt.SNP$same.ALT == 1), ]

  # Compute freq male and freq female
  data.filt.SNP.same.alt$freq.alt.female <- data.filt.SNP.same.alt$count.alt.female/(data.filt.SNP.same.alt$count.ref.female+data.filt.SNP.same.alt$count.alt.female)
  data.filt.SNP.same.alt$freq.alt.male <- data.filt.SNP.same.alt$count.alt.male/(data.filt.SNP.same.alt$count.ref.male+data.filt.SNP.same.alt$count.alt.male)

  ##################
  ##### Gstats #####
  ##################

  # Compute Gstat
  pseudo.raw <- data.filt.SNP.same.alt[, c("count.ref.male", "count.alt.male", "count.ref.female", "count.alt.female")]
  pseudo.raw <- pseudo.raw+1
  gstat.pseudo.raw <- ComputeGstat(input.df = pseudo.raw) 
  data.filt.SNP.same.alt$gstat.pseudo.raw <- gstat.pseudo.raw 

  # Compute smooth Gstat per contig
  w <- 1000
  uni.ctg <- unique(data.filt.SNP.same.alt$CHROM)
  data.filt.SNP.same.alt$gsmooth.pseudo.raw <- NA
  for(ctg in uni.ctg) {
    ind <- which(data.filt.SNP.same.alt$CHROM == ctg)
    tmp <- data.filt.SNP.same.alt[ind, c("POS", "gstat.pseudo.raw")]
    gsmooth <- ksmooth(tmp$POS, tmp$gstat.pseudo.raw, kernel = "normal", bandwidth = w, range.x = range(tmp$POS), n.points = length(tmp$POS), x.points = tmp$POS)
    data.filt.SNP.same.alt$gsmooth.pseudo.raw[ind] <- gsmooth$y
  }

  # Compute H0 
  WGsmooth <- log(data.filt.SNP.same.alt$gsmooth.pseudo.raw)
  median.WGsmooth <- median(WGsmooth)
  sW <- median(abs(WGsmooth[which(WGsmooth<=median.WGsmooth)] - median.WGsmooth))
  ind.outlier <- which((WGsmooth-median.WGsmooth)>(5.2*sW))
  mu <- log(median(data.filt.SNP.same.alt$gsmooth.pseudo.raw[-ind.outlier]))
  mode <- mlv(data.filt.SNP.same.alt$gsmooth.pseudo.raw[-ind.outlier])
  sd <- sqrt(mu - log(mode$M))

  # Test to simulate H0 values 
  pdf("results/call_var/2017_11_14_Gprime_statistics_obs_simulated_underH0.pdf")
  hist(rlnorm(dim(data.filt.SNP.same.alt)[1], mu, sd), breaks = 50, xlim = c(0, 10), main = "Histogram of smoothed G' statistics", xlab = "G'")
  hist(data.filt.SNP.same.alt$gsmooth.pseudo.raw, breaks = 500, add = T, col = adjustcolor("blue", 0.2))
  legend("topright", c(paste("simulated G' statistics according to logN(mu=", round(mu,1), ", sd=", round(sd,1), ")", sep = ""), paste("observed G' statistics (max=", round(max(data.filt.SNP.same.alt$gsmooth.pseudo.raw), 1), " - median=", round(median(data.filt.SNP.same.alt$gsmooth.pseudo.raw), 1), ")", sep = "")), fill = c("white", "blue"), cex = 0.8, bty = "n")
  dev.off()

  # Compute z score
  GZ <- (log(data.filt.SNP.same.alt$gsmooth.pseudo.raw) - mu)/sd
  data.filt.SNP.same.alt$Zscore.gsmooth.pseudo.raw <- GZ

  # Compute pvalue
  data.filt.SNP.same.alt$pval.Zscore.gsmooth.pseudo.raw <- 2*pnorm(-abs(GZ))

  # Adjust for multiple testing
  qobj <- qvalue(p = data.filt.SNP.same.alt$pval.Zscore.gsmooth.pseudo.raw)
  data.filt.SNP.same.alt$qval.Zscore.gsmooth.pseudo.raw <- qobj$qvalues
  pdf("results/call_var/2017_11_14_qstatistics_pvalue_Gprime_statistics_obs.pdf")
  h <- hist(qobj)
  print(h)
  dev.off()
 
  ##################
  ##### Fisher #####
  ##################

  # Add raw pvalue associated to Fisher test done on raw counts
  pval.fisher <- ComputePvalueFisher(input.df = data.filt.SNP.same.alt[, c("count.ref.male", "count.alt.male", "count.ref.female", "count.alt.female")])
  data.filt.SNP.same.alt$pval.fisher.raw.count <- pval.fisher

  # Adjust Fisher pvalue with BH
  data.filt.SNP.same.alt$pval.adjBH.fisher.raw.count <- p.adjust(data.filt.SNP.same.alt$pval.fisher.raw.count, method = "BH")

  write.table(data.filt.SNP.same.alt, "results/call_var/2017_11_14_merge.sexe.common.biallelic.filt.SNP.same.alt.fisher.G.txt", sep = "\t",  quote = F, row.names = F)
}

do.H0 <- F
##### Compute H0 distribution of count tables: UPDATE: instead using Bastide et al Methods to choose the best alpha value, I use instead a vglm fit on betabinomial to get alpha and beta parameters.
tests <- read.csv("results/call_var/2017_11_14_merge.sexe.common.biallelic.filt.SNP.same.alt.fisher.FDR.G.txt", sep = "\t", h = T, stringsAsFactors = F)
fit <- vglm(cbind(c(tests$count.ref.female, tests$count.alt.female), c(tests$count.ref.male, tests$count.alt.male)) ~ 1, betabinomial, trace = T)
coeff <- Coef(fit)
alpha <- (1-coeff["rho"])/(coeff["rho"]*(1+(1-coeff["mu"])/coeff["mu"]))
beta <- alpha*(1-coeff["mu"])/coeff["mu"]


##### Compute 10 datasets under H0 with alpha set to obtain above
do.simu.10x <- F
if(do.simu.10x) {
	data.tests <- read.csv("results/call_var/2017_11_14_merge.sexe.common.biallelic.filt.SNP.same.alt.fisher.G.txt", sep = "\t", h = T, stringsAsFactors = F)
	p.male <- mean(data.tests$tot.male)/(mean(data.tests$tot.female)+mean(data.tests$tot.male))
	alpha.choose <- alpha
	beta.choose <- beta
	ind.ord <- order(data.tests$pval.fisher.raw.count)
	pval.obs.ord <- data.tests$pval.fisher.raw.count[ind.ord]
	data.tests.simu <- data.tests[ind.ord, ]
	for (i in 1:10) {
		print(i)
		cores=detectCores()
		cl <- makeCluster(cores[1]-1) #not to overload your computer
		registerDoParallel(cl)
		pvalH0_out <- foreach(i=1:dim(data.tests)[1], .combine=rbind) %dopar% {
	          	tmp <- data.tests[i, ]
	
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
	        pval.fisher.H0 <- ComputePvalueFisher(input.df = pvalH0_out[, c("count.ref.male", "count.alt.male", "count.ref.female", "count.alt.female")])
	        pvalH0_out$pval.fisher.H0 <- pval.fisher.H0

		saveRDS(pvalH0_out, paste("results/fisher/2017_11_29_simulated_#", i, "_count_H0_alpha", alpha.choose, ".RData", sep = ""))	
	}

	pH0.simu <- NULL
	for (i in 1:10) {
		pvalH0_out <- readRDS(paste("results/fisher/2017_11_15_simulated_#", i, "_count_H0_alpha", alpha.choose, ".RData", sep = ""))$pval.fisher.H0
		pH0.simu <- cbind(pH0.simu, pvalH0_out)
	}
	pH0.simu <- cbind(pH0.simu, data.frame(pval.obs.ord = pval.obs.ord))
	# For the ith rank SNP, # nb simu SNP with pvalue <= pvalue obs for Fisher and compute fdr with rank
	pemp <- apply(pH0.simu, 1, function(x) {y <- unlist(x); length(which(y[1:10] <= y[11]))})/10
	fdr <- p.adjust(pemp, method = "BH")
	data.tests.simu$FDR <- fdr
	write.table(data.tests.simu, "results/call_var/2017_11_14_merge.sexe.common.biallelic.filt.SNP.same.alt.fisher.FDR.G.txt", sep = "\t", col.names = T, row.names = F)
}	

##### Zoom on tests output
size.genome.mbelari <- fread("data/ReferenceGenomes/2017_09_13_Mbelari.sizes.genome", sep = "\t", h = F, stringsAsFactors = F)
genes.pos <- read.csv("data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_genes.bed", h = F, sep = "\t", stringsAsFactors = F)
all <- NULL
tests$gene <- ""
for (ctg in unique(tests$CHROM)) {
	tmp.pos <- genes.pos[which(genes.pos$V1 == ctg), ]
	tmp.tests <- tests[which(tests$CHROM == ctg), ]
	for (i in 1:dim(tmp.tests)[1]){
		pos <- tmp.tests[i, ]$POS
		g = tmp.pos[which(tmp.pos$V2<= pos & pos <= tmp.pos$V3), ]$V4
		if(length(g)>0){
			tmp.tests[i, ]$gene <- paste(g, collapse = "-")
		}
		all <- rbind(all, tmp.tests[i, ])
	}
}

tab <- table(tests$qval.Zscore.gsmooth.pseudo.raw<0.05, tests$CHROM)
tab.signif <- tab[, which(tab[2, ]>0)]
pdf("results/call_var/2017_11_26_signif_G_prime.pdf")
barplot(tab.signif[2, ], xlab = "At least 1 SNP signif. on this contig (qval of G')", ylab = "#signif per contig", xaxt = "n",
main = paste("Output of signif. G' tests (", sum(tab.signif[2, ]), ") per\nSNP at a threshold of 5% per contig (", dim(tab.signif)[2], ")", sep = ""))
dev.off()

tab <- table(tests$FDR<0.05, tests$CHROM)
tab.signif <- tab[, which(tab[2, ]>0)]
pdf("results/call_var/2017_11_26_signif_Fisher.pdf")
barplot(tab.signif[2, ], xlab = "At least 1 SNP signif. on this contig (FDR for fisher)", ylab = "#signif per contig", xaxt = "n",
main = paste("Output of signif. fisher tests (", sum(tab.signif[2, ]), ") per\nSNP at a threshold of 5% per contig (", dim(tab.signif)[2], ")", sep = ""))
dev.off()


##### Plot FDR representation by contig with highest nb of significant SNPs to lowest number of significant SNPs
fisher <- tests[, c("CHROM", "POS", "POS", "REF", "pval.fisher.raw.count", "pval.adjBH.fisher.raw.count", "FDR")]
tab <- table(tests$FDR<0.05, tests$CHROM)
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
