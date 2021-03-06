# Estimate fragment length distribution on processed reads.  
# done on local computer

require(sleuth)

kall <- read_kallisto_h5("results/kallisto_fragmentlength/2017_08_25_abundance.h5")
flen <- unlist(sapply(1:1000, function(x) rep(x, kall$fld[x])))
pdf("results/kallisto_fragmentlength/2017_08_25_hist_fragments_length.pdf")
hist(flen, main = "Histogram of estimated fragment length from assembly only", sub = paste(paste(names(summary(flen)), 
  ":", summary(flen), sep = ""), collapse = " "))
dev.off()
