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
ComputeNormalizedCount <- function(counts.female, counts.male, level = "gene"){
  if(level == "bp") {
    counts.genes <- merge(counts.female[, paste("V", c(8, 1, 6,7,9,4), sep = "")], counts.male[, paste("V", c(4, 11), sep = "")], by = "V8")
    colnames(counts.genes)[6] <- "counts.raw.female"; colnames(counts.genes)[7] <- "counts.raw.male"
  }
  if(level == "gene"){
    counts.genes <- merge(counts.female[, paste("V", c(1:4, 6, 11), sep = "")], counts.male[, paste("V", c(4, 11), sep = "")], by = "V4")
    colnames(counts.genes)[6] <- "counts.raw.female"; colnames(counts.genes)[7] <- "counts.raw.male"
  }
  if(level == "contig"){
    counts.genes <- merge(counts.female[, c("Contig", "Length", "Reads")], counts.male[, c("Contig", "Reads")], by = "Contig")
    colnames(counts.genes)[3] <- "counts.raw.female"; colnames(counts.genes)[4] <- "counts.raw.male"
  }
  counts.genes.to.norm <- matrix(c(counts.genes$counts.raw.female, counts.genes$counts.raw.male), byrow = F, ncol = 2) 
  counts.genes.norm <- normalize.quantiles(counts.genes.to.norm)
  counts.genes <- as.data.frame(counts.genes)
  counts.genes$counts.norm.female <- counts.genes.norm[, 1]
  counts.genes$counts.norm.male <- counts.genes.norm[, 2]
  counts.genes$log2.raw.FC <- log2((counts.genes$counts.raw.female+1)/(counts.genes$counts.raw.male+1))
  counts.genes$log2.norm.FC <- log2((counts.genes$counts.norm.female+1)/(counts.genes$counts.norm.male+1))
  return(counts.genes)
}

##### ComputePvalueFisher #####
ComputePvalueFisher<- function(input.df = merge.sexe.filt[, c("count.ref.male", "count.alt.male", "count.ref.female", "count.alt.female")]) {
  pval.fisher <- apply(input.df, 1, function(x) {
    f = fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol = 2, byrow = T));
    pval = f$p.value;
  })
  return(pval.fisher)
}

##### ComputeGstat #####
ComputeGstat<- function(input.df = merge.sexe.filt[, c("count.ref.male", "count.alt.male", "count.ref.female", "count.alt.female")]) {
  G.stat <- apply(input.df, 1, function(x) {
    tab = matrix(c(x[1], x[2], x[3], x[4]), ncol = 2, byrow = T);
    margin.row1 = x[1] + x[2];
    margin.row2 = x[3] + x[4];
    margin.col1 = x[1] + x[3];
    margin.col2 = x[2] + x[4];
    tot = margin.row1 + margin.row2;
    tab.pred = matrix(c(margin.row1*margin.col1/tot, margin.row1*margin.col2/tot, margin.row2*margin.col1/tot, margin.row2*margin.col2/tot), ncol = 2, byrow = T);
    G = 2*(tab[1,1]*log(tab[1,1]/tab.pred[1,1]) + tab[1,2]*log(tab[1,2]/tab.pred[1,2]) + tab[2,1]*log(tab[2,1]/tab.pred[2,1]) + tab[2,2]*log(tab[2,2]/tab.pred[2,2]));
  })
  return(G.stat)
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
