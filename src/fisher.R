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

