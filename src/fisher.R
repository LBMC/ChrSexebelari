library(gridExtra)
require(modeest)
require(qvalue)
require(foreach)
require(doParallel)
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
    #plot(tmp$POS, tmp$gstat.pseudo.raw)
    #window <- c(500, 1000, 5000, 10000)
    #for (j in seq_along(window)) {
    #  w <- window[j]
    #  pred <- ksmooth(tmp$POS, tmp$gstat.pseudo.raw, kernel = "normal", bandwidth = w, range.x = range(tmp$POS), n.points = length(tmp$POS), x.points = tmp$POS)
    #  points(pred$x, pred$y, type = "l", col = j, lty = 2)
    #}
    #legend("topright", legend = window, fill = seq_along(window))
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
 
  # Add raw pvalue associated to Fisher test done on raw counts
  pval.fisher <- ComputePvalueFisher(input.df = data.filt.SNP.same.alt[, c("count.ref.male", "count.alt.male", "count.ref.female", "count.alt.female")])
  data.filt.SNP.same.alt$pval.fisher.raw.count <- pval.fisher

  # Adjust Fisher pvalue with BH
  data.filt.SNP.same.alt$pval.adjBH.fisher.raw.count <- p.adjust(data.filt.SNP.same.alt$pval.fisher.raw.count, method = "BH")

  write.table(data.filt.SNP.same.alt, "results/call_var/2017_11_14_merge.sexe.common.biallelic.filt.SNP.same.alt.fisher.G.txt", sep = "\t",  quote = F, row.names = F)
}


##### Compute H0 distrib of count tables from Bastide et al, 2013. alpha != beta, alpha tested from 10 to 40

do.H0 <- F
if(do.H0){
	data.tests <- read.csv("results/call_var/2017_11_14_merge.sexe.common.biallelic.filt.SNP.same.alt.fisher.G.txt", sep = "\t", h = T, stringsAsFactors = F)
	p.male <- mean(data.tests$tot.male)/(mean(data.tests$tot.female)+mean(data.tests$tot.male))
	alpha <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 55, 60, 70, 75, 80)
	for (alphanull in alpha) {
		cores=detectCores()
		cl <- makeCluster(cores[1]-1) #not to overload your computer
		registerDoParallel(cl)
		pvalH0_out <- foreach(i=1:dim(data.tests)[1], .combine=rbind) %dopar% {
	          	tmp <- data.tests[i, ]
	
			nref <- tmp$count.ref.female + tmp$count.ref.male
	                betanull <- (alphanull-alphanull*p.male)/p.male
	                pmref <- rbeta(1, shape1 = alphanull, shape2 = betanull)
			nrefmale <- rbinom(n = 1, size = nref, prob = pmref)
			nreffemale <- nref-nrefmale
	
			nalt <- tmp$count.alt.female + tmp$count.alt.male
                	pmalt <- rbeta(1, shape1 = alphanull, shape2 = betanull)
			naltmale <- rbinom(n = 1, size = nalt, prob = pmalt)
			naltfemale <- nalt-naltmale

 			data.frame(ID = tmp$ID, count.ref.male = nrefmale, count.alt.male = naltmale, count.ref.female = nreffemale, count.alt.female = naltfemale, stringsAsFactors = F)
		}
	  	stopCluster(cl)
	        pval.fisher.H0 <- ComputePvalueFisher(input.df = pvalH0_out[, c("count.ref.male", "count.alt.male", "count.ref.female", "count.alt.female")])
	        pvalH0_out$pval.fisher.H0 <- pval.fisher.H0

		saveRDS(pvalH0_out, paste("results/fisher/2017_11_14_simulated_count_H0_from_Bastide_alpha", alphanull, "_biallelic.filt.SNP.same.alt.fisher.RData", sep = ""))	
	}

	# Choose alpha
        h <- hist(data.tests$pval.fisher.raw.count, breaks = seq(0, 1, length = 700), plot = F, include.lowest = T, right = T)
	p.bins <- h$counts/dim(data.tests)[1]
	res.choose.alpha <- NULL
	for (alphanull in alpha) {
		tmp <- readRDS(paste("results/fisher/2017_11_14_simulated_count_H0_from_Bastide_alpha", alphanull, "_biallelic.filt.SNP.same.alt.fisher.RData", sep = ""))

        	htmp <- hist(tmp$pval.fisher.H0, seq(0, 1, length = 700), plot = F, include.lowest = T, right = T)
		ptest.chisq <- chisq.test(htmp$counts, p = p.bins, correct = F)$p.value
		ptest.fisher <- fisher.test(cbind(htmp$counts, h$counts), simulate.p.value = T)$p.value
		ptest.ks <- ks.test(htmp$counts, h$counts)$p.value

		# Quantile-quantile plot of observed vs. simulated P-values for the Fisher’s exact test on allele frequency differentiation
		png(paste("results/fisher/2017_11_14_qqplot_simulated_vs_obs_pvalues_alpha", alphanull, ".png", sep = ""), w = 850)
		par(mfrow = c(1,2))
		qqplot(tmp$pval.fisher.H0, data.tests$pval.fisher.raw.count, main = paste("Quantile-quantile plot of observed vs. simulated\nP-values for the Fisher exact test - alpha=", alphanull, " on 700 bins", sep = ""), xlab = "Simulated null pvalue", ylab = "Observed pvalue")
		text(0.6, 0.4, paste("pchisq=", ptest.chisq))
		text(0.6, 0.2, paste("pfisher=", ptest.fisher))
		text(0.6, 0.0, paste("pks=", ptest.ks))
		abline(a=0,b=1,lty=2)
		qqplot(-log10(tmp$pval.fisher.H0), -log10(data.tests$pval.fisher.raw.count), main = paste("Quantile-quantile plot of observed vs. simulated\n-log10(P-values) for the Fisher exact test - alpha=", alphanull, sep = ""), xlab = "Simulated null -log10(pvalue)", ylab = "Observed -log10(pvalue)")
		abline(a=0,b=1,lty=2)
		dev.off()

		res.choose.alpha <- rbind(res.choose.alpha, data.frame(alpha = alphanull, beta = (alphanull-alphanull*p.male)/p.male, pchisq_pvalH0_to_obs = ptest.chisq, pfisher_pvalH0_to_obs = ptest.fisher, pks_on_counts_pvalH0_to_obs = ptest.ks, stringsAsFactors = F))
	}
	write.table(res.choose.alpha, "results/fisher/2017_11_14_summary_choose_alpha_for_siimulated_count_H0_from_Bastide_biallelic.filt.SNP.same.alt.fisher.txt", sep = "\t", col.names = T, row.names = F, quote = F)
}
# Maybe test on 1000bins as well

##### Compute 10 datasets under H0 with alpha = 75
do.simu.10x <- F
if(do.simu.10x) {
	data.tests <- read.csv("results/call_var/2017_11_14_merge.sexe.common.biallelic.filt.SNP.same.alt.fisher.G.txt", sep = "\t", h = T, stringsAsFactors = F)
	p.male <- mean(data.tests$tot.male)/(mean(data.tests$tot.female)+mean(data.tests$tot.male))
	alpha.choose <- 75
	beta.choose <- (alpha.choose-alpha.choose*p.male)/p.male
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
	                beta.choose <- (alpha.choose-alpha.choose*p.male)/p.male
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

		saveRDS(pvalH0_out, paste("results/fisher/2017_11_15_simulated_#", i, "_count_H0_alpha", alpha.choose, ".RData", sep = ""))	

		# For the ith rank SNP, # nb simu SNP with pvalue < pvalue obs for Fisher
		tmp <- data.frame(nb.simu.SNPs.lower.pval = sapply(pval.obs.ord, function(x) length(which(pvalH0_out$pval.fisher.H0<x))), stringsAsFactors = F)
		saveRDS(tmp, paste("results/fisher/2017_11_15_simulated_compare_to_obs_#", i, "_count_H0_alpha", alpha.choose, ".RData", sep = ""))	
	}
}



