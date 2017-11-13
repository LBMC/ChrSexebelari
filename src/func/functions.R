##### ExtractCountFromVCF #####
ExtractCountFromVCF <- function(vcf.without.header, sexe){
  info <- unname(do.call(rbind, sapply(vcf.without.header$INFO, function(x) {
    eval(parse(text = paste("DP4=c(", strsplit(strsplit(x, "DP4=")[[1]][2], ";")[[1]][1], ")", sep = "")));
    is.INDEL <- length(grep("INDEL", x));
    ref = sum(DP4[1:2]);
    alt = sum(DP4[3:4]);
    tot=sum(DP4);
    AN = strsplit(strsplit(x, "AN=")[[1]][2], ";")[[1]][1];
    as.numeric(c(AN, ref, alt, tot, is.INDEL))
  }, simplify = F)))
  info <- as.data.frame(info)
  colnames(info) <- c("nb.alleles", paste("count.ref.", sexe, sep = ""), paste("count.alt.", sexe, sep = ""), paste("tot.", sexe, sep = ""), "is.INDEL")
  return(info)
}

##### ComputeNormalizedCount #####
ComputeNormalizedCount <- function() {

}

##### ComputePvalueFisher #####
ComputePvalueFisher<- function(input.df = merge.sexe.filt[, c("count.ref.male", "count.alt.male", "count.ref.female", "count.alt.female")]) {
  pval.fisher <- apply(input.df, 1, function(x) {
    f = fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol = 2, byrow = T));
    pval = f$p.value;
  })
  return(pval.fisher)
}


##### SummarizeSNPsINDELsWithinDataFrame #####
SummarizeSNPsINDELsWithinDataFrame <- function(merge.df, suff = "shared position") {
  ind.genic <- which(merge.df$genes != "")
  df <- data.frame(condition = suff, nb.SNPs.all = table(merge.df$is.INDEL)["0"],
  nb.INDELs.all =  table(merge.df$is.INDEL)["1"],
  nb.SNPs.genic = table(merge.df[ind.genic]$is.INDEL)["0"],
  nb.INDELs.genic = table(merge.df[ind.genic]$is.INDEL)["1"],
  nb.ctg.SNPs.all = length(unique(merge.df$CHROM[which(merge.df$is.INDEL == 0)])),
  nb.ctg.INDELs.all = length(unique(merge.df$CHROM[which(merge.df$is.INDEL == 1)])), stringsAsFactors = F)
  return(df)
}
