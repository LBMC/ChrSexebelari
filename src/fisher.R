require(data.table)
require(dplyr)
library(foreach)
library(doParallel)
library(gridExtra)

do.fisher <- F
alpha <- 0.05

if (do.fisher) {
  data.filt.SNP <- read.csv("results/call_var/2017_11_09_merge.sexe.common.biallelic.filt.SNP.txt", sep = "\t", h = T, stringsAsFactors = F)
  # Compute normalized counts by median number of counts and 1000 individuals at a given position (DP field) 
  median.tot.male <- median(data.filt.SNP$tot.male); median.tot.female <- median(data.filt.SNP$tot.female);

  #norm.male <- ComputeNormalizedCount(counts.ref = data.filt.SNP$count.ref.male, counts.alt = data.filt.SNP$count.ref.male, norm.factor = median.tot.male) 
  #data.filt.SNP$count.norm.ref.male <- norm.male[["norm.ref"]]
  #data.filt.SNP$count.norm.alt.male <- norm.male[["norm.alt"]]

  #norm.female <- ComputeNormalizedCount(counts.ref = data.filt.SNP$count.ref.female, counts.alt = data.filt.SNP$count.ref.female, norm.factor = median.tot.female) 
  #data.filt.SNP$count.norm.ref.female <- norm.female[["norm.ref"]]
  #data.filt.SNP$count.norm.alt.female <- norm.female[["norm.alt"]]
  
  # Add raw pvalue associated to Fisher test done on raw counts
  pval.fisher <- ComputePvalueFisher(input.df = data.filt.SNP[, c("count.ref.male", "count.alt.male", "count.ref.female", "count.alt.female")])
  data.filt.SNP$pval.fisher.raw.count <- pval.fisher
  
  #pval.fisher.norm <- ComputePvalueFisher(input.df = data.filt.SNP[, c("count.orm.ref.male", "count.norm.alt.male", "count.norm.ref.female", "count.norm.alt.female")])
  #data.filt.SNP$pval.fisher.norm.count <- pval.fisher.norm

  # Compute freq male and freq female
  data.filt.SNP$freq.alt.female <- data.filt.SNP$count.alt.female/(data.filt.SNP$count.ref.female+data.filt.SNP$count.alt.female)
  data.filt.SNP$freq.alt.male <- data.filt.SNP$count.alt.male/(data.filt.SNP$count.ref.male+data.filt.SNP$count.alt.male)
  #data.filt.SNP$freq.norm.alt.female <- data.filt.SNP$count.norm.alt.female/(data.filt.SNP$count.norm.ref.female+data.filt.SNP$count.norm.alt.female)
  #data.filt.SNP$freq.norm.alt.male <- data.filt.SNP$count.norm.alt.male/(data.filt.SNP$count.norm.ref.male+data.filt.SNP$count.norm.alt.male)

  # Adjust pvalue with BH
  data.filt.SNP$pval.adjBH.fisher.raw.count <- p.adjust(data.filt.SNP$pval.fisher.raw.count, method = "BH")
  #data.filt.SNP$pval.adjBH.fisher.norm.count <- p.adjust(data.filt.SNP$pval.fisher.norm.count, method = "BH")

  write.table(data.filt.SNP, "results/call_var/2017_11_10_merge.sexe.common.biallelic.filt.SNP.fisher_raw_count.txt", sep = "\t",  quote = F, row.names = F)
}
