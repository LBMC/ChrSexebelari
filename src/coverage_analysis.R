require(data.table)
require(dplyr)
require(ggplot2)
library(devtools)
library(Biobase)
library(preprocessCore)
require(MASS)
require(foreach)
require(doParallel)
require(parallel)
source("src/func/functions.R")
  
# Generate gff3 containing gene only
gff.mbela <- read.csv("data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2.gff3", 
  sep = "\t", h = F, stringsAsFactors = F)
gff.genes.mbela <- gff.mbela[which(gff.mbela$V3 == "gene"), ]
write.table(gff.genes.mbela, "data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_genes.gff3", 
  col.names = F, row.names = F, quote = F, sep = "\t")

contigs.mbela <- fread("data/ReferenceGenomes/2017_09_13_Mbelari.sizes.genome", sep = "\t", 
  h = F, stringsAsFactors = F)

thresh.norm.FC <- 1

compute.norm <- F

if(compute.norm){
  ##### Plot cluster size when comparing unmapped female and male reads
  clus.size <- read.csv("results/mapping/unmapped/2017_11_22_cdhit.fa.clstr.size", 
    sep = "\t", h = T, stringsAsFactors = F)
  # > dim(clus.size) [1] 9079667 2 summary(clus.size$nb) Min.  1st Qu.  Median Mean
  # 3rd Qu.  Max.  0.000 0.000 0.000 0.713 0.000 4430.000
  pdf("results/mapping/unmapped/cluster_size.pdf")
  hist(log10(clus.size$nb), breaks = 50, main = paste("log10(size of clusters) (n=", 
    dim(clus.size)[1], ") obtained with CD-HIT\nto compare unmapped female and male reads", 
    sep = ""), cex.main = 0.9)
  dev.off()
  system("bash src/date.sh results/mapping/mapped/cluster_size.pdf")

  ##### Pick contig name and length
  print(summary(contigs.mbela$V2))

  #####  Percentage of genic region:
  tot.size <- sum(contigs.mbela$V2)
  gene.size <- abs(gff.genes.mbela$V5 - gff.genes.mbela$V4)
  perc.genic <- 100 * sum(gene.size)/tot.size  # 52.72679
  perc.ctg <- sapply(1:dim(contigs.mbela)[1], function(x) {
    ctg = contigs.mbela[x]
    tmp = gff.genes.mbela[which(gff.genes.mbela$V1 == ctg), ]
    s = sum(abs(tmp$V5 - tmp$V4))
    100 * s/contigs.mbela$V2[x]
  })
  # > summary(perc.ctg) Min.  1st Qu.  Median Mean 3rd Qu.  Max.  0.000000 0.000000
  # 0.000000 0.000373 0.000000 6.082403

  #####  Summary nb genes per contig
  gff.genes.mbela$V10 <- sapply(gff.genes.mbela$V9, function(x) strsplit(x, "ID=")[[1]][2])
  genes.contigs.mbela <- table(gff.mbela$V1, gff.mbela$V3)
  pdf("results/coverage/number_genes_detected_in_contigs.pdf")
  barplot(table(genes.contigs.mbela[, "gene"]), main = paste("nb genes within contig (genic% = ", 
    perc.genic, ")", sep = ""), ylab = "Frequency", xlab = "nb genes per contig")
  dev.off()

  ##### Nb of reads per contig
  cov.female <- fread("results/coverage/2017_10_26_coverage_summary_tablet_female_JU2817_trim_Mbelari_mapped_rmdup_rg_realign_indels.txt", 
    sep = "\t", h = T, stringsAsFactors = F)
  cov.male <- fread("results/coverage/2017_10_26_coverage_summary_tablet_male_JU2817_trim_Mbelari_mapped_rmdup_rg_realign_indels.txt", 
    sep = "\t", h = T, stringsAsFactors = F)

  ##### Summary number of reads on contigs per pool after having removed duplicates
  cov.female <- cov.female[which(cov.female$Contig %in% contigs.mbela$V1), ]
  cov.male <- cov.male[which(cov.male$Contig %in% contigs.mbela$V1), ]
  print(summary(cov.female$Reads))
  # Min. 1st Qu.  Median Mean 3rd Qu.  Max.  0 108 478 2170 2455 108520
  print(summary(cov.male$Reads))
  # Min. 1st Qu.  Median Mean 3rd Qu.  Max.  0 62 284 1230 1351 101118

  ##### Mean coverage value per sample
  print(sum(cov.male$Reads * 100 * 2)/sum(cov.male$Length))  #24.76252
  print(sum(cov.female$Reads * 100 * 2)/sum(cov.female$Length))  #43.67202

  ##### Normalize genes counts  Generate bed file for gene region coordinates
  system("bash src/convertgff_to_bed.sh data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_genes.gff3 data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_genes.bed")
  source("src/format_bed.R")
  # At gene bp level
  system("bash src/intersect_bed_genes.sh")
  # At gene level
  system("bash src/get_counts_per_gene.sh")

  ##### Normalize counts at contig level, gene level and at gene bp level
  # at contig
  counts.contigs <- ComputeNormalizedCount(cov.female, cov.male, "contig")
  counts.contigs$missing.sexe <- ""; counts.contigs$missing.sexe[which(counts.contigs$counts.raw.female == 0)] <- "female"; 
  counts.contigs$missing.sexe[which(counts.contigs$counts.raw.male == 0)] <- "male"; 
  write.table(counts.contigs, "results/coverage/FC_normalized_coverage_at_contig.txt", 
    sep = "\t", quote = F, col.names = T, row.names = F)
  system("bash src/date.sh results/coverage/FC_normalized_coverage_at_contig.txt")

  # at gene 
  counts.genes.female <- fread("results/coverage/2017_12_06_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_count_genes.txt", 
    h = F, sep = "\t", stringsAsFactors = F)
  counts.genes.male <- fread("results/coverage/2017_12_06_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_count_genes.txt", 
    h = F, sep = "\t", stringsAsFactors = F)
  counts.genes <- ComputeNormalizedCount(counts.genes.female, counts.genes.male, "gene")
  counts.genes$missing.sexe <- ""; counts.genes$missing.sexe[which(counts.genes$counts.raw.female == 0)] <- "female"; 
  counts.genes$missing.sexe[which(counts.genes$counts.raw.male == 0)] <- "male"; 
  write.table(counts.genes, "results/coverage/FC_normalized_coverage_at_gene.txt", 
    sep = "\t", quote = F, col.names = T, row.names = F)
  system("bash src/date.sh results/coverage/FC_normalized_coverage_at_gene.txt")

  ##### Intersect genes and coverage at bp resolution
  counts.genes.female <- tbl_df(fread("results/coverage/2017_12_06_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_sort_in_genes.bed", 
    h = F, sep = "\t", stringsAsFactors = F))
  counts.genes.male <- tbl_df(fread("results/coverage/2017_12_06_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_sort_in_genes.bed", 
    h = F, sep = "\t", stringsAsFactors = F))
  # in_male <- setdiff(unique(counts.genes.male$V8),
  # unique(counts.genes.female$V8)) in_female <-
  # setdiff(unique(counts.genes.female$V8), unique(counts.genes.male$V8)) >
  # in_female #24211 genes [1] 'MBELA.g12439' 'MBELA.g12585' 'MBELA.g14323'
  # 'MBELA.g14395' 'MBELA.g14542' [6] 'MBELA.g17943' 'MBELA.g23060' > in_male
  # #24209 genes [1] 'MBELA.g11678' 'MBELA.g15899' 'MBELA.g17246' 'MBELA.g17397'
  # 'MBELA.g20649'

  ##### Compute subset of data to do merge
  counts.genes.bp <- ComputeNormalizedCount(counts.genes.female, counts.genes.male, 
    "bp")
  write.table(counts.genes.bp, "results/coverage/FC_normalized_coverage_at_bp_within_genes.txt", 
    sep = "\t", quote = F, col.names = T, row.names = F)
  system("bash src/date.sh results/coverage/FC_normalized_coverage_at_bp_within_genes.txt")
}

##### Zoom on gene with no count in one gene or in one contig
counts.contigs <- read.csv("results/coverage/2017_12_06_FC_normalized_coverage_at_contig.txt", 
  sep = "\t", h = T, stringsAsFactors = F)
counts.genes <- read.csv("results/coverage/2017_12_06_FC_normalized_coverage_at_gene.txt", 
  sep = "\t", h = T, stringsAsFactors = F)

counts.genes.missing.one.sexe <- counts.genes[which(counts.genes$missing.sexe != ""), ]
counts.genes.missing.one.sexe <- counts.genes.missing.one.sexe[order(counts.genes.missing.one.sexe$missing.sexe), ]
counts.contigs.missing.one.sexe <- counts.contigs[which(counts.contigs$missing.sexe != ""), ]
counts.contigs.missing.one.sexe <- counts.contigs.missing.one.sexe[order(counts.contigs.missing.one.sexe$missing.sexe), ]
write.table(counts.genes.missing.one.sexe, "results/coverage/counts_per_genes_raw_counts_not_present_in_both_sexe.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(counts.contigs.missing.one.sexe, "results/coverage/counts_per_contigs_raw_counts_not_present_in_both_sexe.txt", sep = "\t", col.names = T, row.names = F, quote = F)
system("bash src/date.sh results/coverage/counts_per_contigs_raw_counts_not_present_in_both_sexe.txt")
system("bash src/date.sh results/coverage/counts_per_genes_raw_counts_not_present_in_both_sexe.txt")

do.test.at.bp <- F
if(do.test.at.bp ){
  ##### In case of bp level, implement test per gene on estimated FC values: if not
  ##### gaussian, compute i) a median like based test on the log2(FC), ii) try the the
  ##### ans comb transformation and perform a test on it.
  counts.bp <- tbl_df(fread("results/coverage/2017_11_30_FC_normalized_coverage_at_bp_within_genes.txt", 
    sep = "\t", stringsAsFactors = F))

  counts.genes.missing.one.sexe <- read.csv("results/coverage/2017_12_21_counts_per_genes_raw_counts_not_present_in_both_sexe.txt", sep = "\t", h = T, stringsAsFactors = F)
  counts.contigs.missing.one.sexe <- read.csv("results/coverage/2017_12_21_counts_per_contigs_raw_counts_not_present_in_both_sexe.txt", sep = "\t", h = T, stringsAsFactors = F)

  #pdf("results/coverage/count_distribution_before_after_cov_at_gene_bp.pdf", w = 12, 
  #  h = 8)
  #par(mfrow = c(2, 2))
  #plot(density(log(counts.bp$counts.raw.male + 1)), main = "Raw male counts at gene bp level", 
  #  xlab = "log+1 raw counts at gene bp level")
  #plot(density(log(counts.bp$counts.raw.female + 1)), main = "Raw female counts at gene bp level", 
  #  xlab = "log+1 raw counts at gene bp level")
  #plot(density(log(counts.bp$counts.norm.male + 1)), main = "Norm. male counts at gene bp level", 
  #  xlab = "log+1 quantile normalized counts at gene bp level")
  #plot(density(log(counts.bp$counts.norm.female + 1)), main = "Norm. female counts at gene bp level", 
  #  xlab = "log+1 quantile normalized counts at gene bp level")
  #dev.off()

  genes <- counts.bp %>% distinct(V8); genes <- unique(genes$V8)
  counts.genes <- read.csv("results/coverage/2017_12_06_FC_normalized_coverage_at_gene.txt", 
    sep = "\t", h = T, stringsAsFactors = F)
  res <- NULL
  for(g in genes) {
    tmp <- as.data.frame(filter(counts.bp, V8 == g)[, c("V1", "counts.norm.female", "counts.norm.male", "log2.norm.FC", "counts.raw.female", "counts.raw.male", "log2.raw.FC")])
    contig <- unique(tmp$V1)
    for (type in c("raw", "norm")) {	
      f <- tmp[, paste("counts.", type, ".female", sep = "")]; f.a <- 2*sqrt(f+3/8)
      m <- tmp[, paste("counts.", type, ".male", sep = "")]; m.a <- 2*sqrt(m+3/8)
      fc <- tmp[, paste("log2.", type, ".FC", sep = "")]; fc.a <- log2(f.a/m.a)
      me <- median(fc); me.a <- median(fc.a)
      mean <- mean(fc); mean.a <- mean(fc.a)
      me.f <- median(f); me.f.a <- median(f.a)
      mean.f <- mean(f); mean.f.a <- mean(f.a)
      me.m <- median(m); me.m.a <- median(m.a) 
      mean.m <- mean(m); mean.m.a <- mean(m.a)
      sd.f <- sd(f); sd.f.a <- sd(f.a)
      sd.m <- sd(m); sd.m.a <- sd(m.a)
      sd.fc <- sd(fc); sd.fc.a <- sd(fc.a)

      if (length(unique(fc)) == 1) {
        med.test <- NA; med.test.a <- NA	
        t.test <- NA; t.test.a <- NA
        t.test.both <- NA; t.test.both.a <- NA
        t.test.both.greater <- NA; t.test.both.a.greater <- NA
        t.test.both.lower <- NA; t.test.both.a.lower <- NA
      } else{
        med.test <- prop.test(sum(fc > 0), length(fc), p = 0.5, "two.sided")$p.value; med.test.a <- prop.test(sum(fc.a > 0), length(fc.a), p = 0.5, "two.sided")$p.value	
        t.test <- t.test(fc)$p.value; t.test.a <- t.test(fc.a)$p.value	
        t.test.both <- t.test(f, m)$p.value; t.test.both.a <- t.test(f.a, m.a)$p.value
        t.test.both.lower <- t.test(f, m, alternative = "less")$p.value; t.test.both.greater <- t.test(f, m, alternative = "greater")$p.value
        t.test.both.a.lower <- t.test(f.a, m.a, alternative = "less")$p.value; t.test.both.a.greater <- t.test(f.a, m.a, alternative = "greater")$p.value
      }
      tmp.res <- data.frame(contig = contig, gene = g, type = type,
        mean.f = mean.f, median.f = me.f, sd.f = sd.f, 
        mean.m = mean.m, median.m = me.m, sd.m = sd.m,
        mean.fc = mean, median.fc = me, sd.fc = sd.fc,
        mean.f.a = mean.f.a, median.f.a = me.f.a, sd.f.a = sd.f.a, 
        mean.m.a = mean.m.a, median.m.a = me.m.a, sd.m.a = sd.m.a,
        mean.fc.a = mean.a, median.fc.a = me.a, sd.fc.a = sd.fc.a,
        pval.med = med.test, pval.med.a = med.test.a, 
        pval.ttest = t.test, pval.ttest.a = t.test.a, 
        pval.ttest.both = t.test.both, pval.ttest.both.a = t.test.both.a, 
        pval.ttest.lower.both = t.test.both.lower, pval.ttest.both.lower.a = t.test.both.a.lower, 
        pval.ttest.greater.both = t.test.both.greater, pval.ttest.both.greater.a = t.test.both.a.greater, 
        stringsAsFactors = F)
      res <- rbind(res, tmp.res)
    }
  }
  res$is.gene.missing.sexe <- sapply(res$gene, function(x) {tmp <- which(counts.genes.missing.one.sexe$V4 == x); ifelse(length(tmp)>0, counts.genes.missing.one.sexe[tmp, "missing.sexe"], "")})
  res$is.contig.missing.sexe <- sapply(res$contig, function(x) {tmp <- which(counts.contigs.missing.one.sexe$Contig == x); ifelse(length(tmp)>0, counts.contigs.missing.one.sexe[tmp, "missing.sexe"], "")})

  write.table(res, "results/coverage/tests_FC_normalized_coverage_at_bp_within_genes.txt", sep= "\t", quote =F, col.names = T, row.names = F)

  system("bash src/date.sh results/coverage/tests_FC_normalized_coverage_at_bp_within_genes.txt")
}

counts.genes.bp <- read.csv("results/coverage/2017_12_20_tests_FC_normalized_coverage_at_bp_within_genes.txt", 
  sep = "\t", h = T, stringsAsFactors = T)
counts.genes.bp.raw <- counts.genes.bp[which(counts.genes.bp$type == "raw"), ]
counts.genes.bp <- counts.genes.bp[which(counts.genes.bp$type == "norm"), ]


##### Investigation on features missing in one sexe
counts.genes.missing.one.sexe <- read.csv("results/coverage/2017_12_21_counts_per_genes_raw_counts_not_present_in_both_sexe.txt", sep = "\t", h = T, stringsAsFactors = F)
counts.contigs.missing.one.sexe <- read.csv("results/coverage/2017_12_21_counts_per_contigs_raw_counts_not_present_in_both_sexe.txt", sep = "\t", h = T, stringsAsFactors = F)

pdf("results/coverage/counts_per_features_raw_counts_not_present_in_both_sexe.pdf")
par(mfrow = c(1, 2))
barplot(table(counts.genes.missing.one.sexe$missing.sexe), main = "#genes missing in one sexe", xlab = "sexe")
barplot(table(counts.contigs.missing.one.sexe$missing.sexe), main = "#contigs missing in one sexe", xlab = "sexe")
dev.off()
system("bash src/date.sh results/coverage/counts_per_features_raw_counts_not_present_in_both_sexe.pdf")

counts.contigs.missing.one.sexe$missing.sexe <- as.factor(counts.contigs.missing.one.sexe$missing.sexe)
counts.genes.missing.one.sexe$missing.sexe <- as.factor(counts.genes.missing.one.sexe$missing.sexe)
graw <- ggplot(counts.genes.missing.one.sexe, aes(x=counts.raw.male, y=counts.raw.female, color=missing.sexe, size=V3-V2)) +
  geom_point(alpha = 0.4) + ggtitle("raw count per gene")+ theme(legend.position="none")
gnorm <- ggplot(counts.genes.missing.one.sexe, aes(x=counts.norm.male, y=counts.norm.female, color=missing.sexe, size=V3-V2)) +
  geom_point(alpha = 0.4) + ggtitle("norm count per gene") + labs(colour = "missin sexe", size = "length (bp)") + 
  theme(legend.justification = "top",legend.text = element_text(size=6))
craw <- ggplot(counts.contigs.missing.one.sexe, aes(x=counts.raw.male, y=counts.raw.female, color=missing.sexe, size=Length)) +
  geom_point(alpha = 0.4) + ggtitle("raw count per contig")+ theme(legend.position="none")
cnorm <- ggplot(counts.contigs.missing.one.sexe, aes(x=counts.norm.male, y=counts.norm.female, color=missing.sexe, size=Length)) +
  geom_point(alpha = 0.4)+ ggtitle("norm count per contig")+ theme(legend.position="none")

pdf("results/coverage/counts_per_features_not_present_in_both_sexe.pdf", w = 9, h = 4)
multiplot(graw, gnorm, craw, cnorm, cols = 2)
dev.off()
system("bash src/date.sh results/coverage/counts_per_features_not_present_in_both_sexe.pdf")

##### Histogram of log2(FC) at contig, gene and estimated at gene bp level
pdf("results/coverage/all_log2_FC.pdf", w = 10, h = 5)
par(mfrow = c(1, 3))
hist(counts.contigs$log2.raw.FC, main = "Histogram of log2(female/male)) for contig", 
  xlab = "log2(FC)", breaks = 100)
hist(counts.contigs$log2.norm.FC, breaks = 100, add = T, col = adjustcolor("red", 
  0.75))
legend("topleft", title = paste("at contig level (n=", dim(counts.contigs)[1], ")", 
  sep = ""), legend = c("log2(FC) on quantile normalized counts", "log2(FC) on raw counts"), 
  fill = c("red", "white"), cex = 0.75, bty = "n")
hist(counts.genes$log2.raw.FC, main = "Histogram of log2(female/male)) for genes", 
  xlab = "log2(FC)", breaks = 100)
hist(counts.genes$log2.norm.FC, breaks = 100, add = T, col = adjustcolor("red", 0.75))
legend("topleft", title = paste("at gene level (n=", dim(counts.genes)[1], ")", sep = ""), 
  legend = c("log2(FC) on quantile normalized counts", "log2(FC) on raw counts"), 
  fill = c("red", "white"), cex = 0.75, bty = "n")
hist(counts.genes.bp.raw$mean.fc, main = "Histogram of log2(female/male)) for\ngenes from bp resolution", 
  xlab = "log2(FC)", breaks = 100)
hist(counts.genes.bp$mean.fc, breaks = 100, add = T, col = adjustcolor("red", 0.75))
legend("topleft", title = paste("at gene bp level (n=", dim(counts.genes)[1], ")", sep = ""), 
  legend = c("log2(FC) on quantile normalized counts", "log2(FC) on raw counts"), 
  fill = c("red", "white"), cex = 0.75, bty = "n")
dev.off()
system("bash src/date.sh results/coverage/all_log2_FC.pdf")

##### Zoom on data with abs(FC)>=threshold 
tab.genes.per.contig <- table(counts.genes$V1)
for (sexe in c("female", "male")) {
  if(sexe == "female"){
    fc.at.bp <- counts.genes.bp[which(counts.genes.bp$mean.fc >= thresh.norm.FC), ]
  }else{
    fc.at.bp <- counts.genes.bp[which(counts.genes.bp$mean.fc <= -thresh.norm.FC), ]  
  }
  tab <- table(fc.at.bp$contig)
  ind <- which(tab <= 1); tab.tmp <- tab[-ind]; lim <- max(as.vector(tab.tmp))
  names(tab.tmp) <- sapply(names(tab.tmp), function(x) paste(x, " (", 
    tab.genes.per.contig[x], ")", sep = ""))
  tab.tmp <- c(tab.tmp, length(which(tab == thresh.norm.FC)))
  names(tab.tmp)[length(tab.tmp)] <- "contigs with1/1 gene"
  tmp.coord <- do.call(rbind, t(sapply(fc.at.bp$gene, function(x) gff.genes.mbela[which(gff.genes.mbela$V9 == paste("ID=", x, sep = "")), c("V4", "V5")], simplify = F)))
  colnames(tmp.coord) <- c("beg", "end")
  info <- cbind(fc.at.bp[, c("contig", "gene")], tmp.coord, fc.at.bp[, c("mean.f", "median.f", "sd.f", "mean.m", "median.m", "sd.m", "mean.fc", "median.fc", "sd.fc")])
  info$contig <- as.character(info$contig)

  eval(parse(text = paste("lim.", sexe, " <- lim", sep = "")))
  eval(parse(text = paste("info.", sexe, " <- info", sep = "")))
  eval(parse(text = paste("fc.at.bp.", sexe, " <- fc.at.bp", sep = "")))
  eval(parse(text = paste("tab.", sexe, ".tmp <- tab.tmp", sep = "")))
}
subset <- rbind(info.female, info.male)
fc.at.bp.male$contig <- as.character(fc.at.bp.male$contig) 
fc.at.bp.female$contig <-  as.character(fc.at.bp.female$contig)
counts.genes.bp$contig <- as.character(counts.genes.bp$contig)

add.others.genes <- do.call(rbind, sapply(unique(c(fc.at.bp.male$contig, fc.at.bp.female$contig)), function(x) counts.genes.bp[which(counts.genes.bp$contig == x & counts.genes.bp$mean.fc > -thresh.norm.FC & counts.genes.bp$mean.fc < thresh.norm.FC), ], simplify = F))
tmp.coord <- do.call(rbind, t(sapply(add.others.genes$gene, function(x) gff.genes.mbela[which(gff.genes.mbela$V9 == paste("ID=", x, sep = "")), c("V4", "V5")], simplify = F)))
colnames(tmp.coord) <- c("beg", "end")
add.others.genes <- cbind(add.others.genes, tmp.coord)
all.subset <- rbind(subset, add.others.genes[, c("contig", "gene", "beg", "end", "mean.f", "median.f", "sd.f", "mean.m", "median.m", "sd.m", "mean.fc", "median.fc", "sd.fc")])

all.subset$contig <- as.factor(all.subset$contig)
all.subset <- all.subset[order(all.subset$contig), ]
all.subset$is.gene.missing.sexe <- sapply(all.subset$gene, function(x) {tmp <- which(counts.genes.missing.one.sexe$V4 == x); 
  ifelse(length(tmp)>0, counts.genes.missing.one.sexe[tmp, "missing.sexe"], "")})
all.subset$is.contig.missing.sexe <- sapply(all.subset$contig, function(x) {tmp <- which(counts.contigs.missing.one.sexe$Contig == x); 
  ifelse(length(tmp)>0, counts.contigs.missing.one.sexe[tmp, "missing.sexe"], "")})
write.table(all.subset, paste("results/coverage/subset_genes_FC_norm_threshold", thresh.norm.FC, ".txt", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)
system(paste("bash src/date.sh results/coverage/subset_genes_FC_norm_threshold", thresh.norm.FC, ".txt", sep = ""))

# barplot 
trans <- function(x){pmin(x,lim) + 0.05*pmax(x-lim,0)}
yticks <- c(round(seq(0, lim.female, length = 3)), seq(lim.female, max(tab.female.tmp)+20, by = 20))
dat <- data.frame(label = names(tab.female.tmp), value = trans(tab.female.tmp))
gf <- ggplot(data=dat, aes(x=label, y=value)) +
  geom_col(position="dodge") +
  geom_rect(aes(xmin=0, xmax=length(tab.female.tmp)+1, ymin=trans(lim.female), ymax=trans(lim.female)+1), fill="white") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(limits=c(0,NA), breaks=trans(yticks), labels=yticks) +
  labs(y="# genes") + ggtitle (paste("#log2(FC)>=", thresh.norm.FC,  " per contig - female enrichment\nContig (nb tot genes)", sep = ""))

yticks <- c(round(seq(0, lim.male, length = 3)), seq(lim.male, max(tab.male.tmp)+20, by = 20))
dat <- data.frame(label = names(tab.male.tmp), value = trans(tab.male.tmp))
gm <- ggplot(data=dat, aes(x=label, y=value)) +
  geom_col(position="dodge") +
  geom_rect(aes(xmin=0, xmax=length(tab.male.tmp)+1, ymin=trans(lim.male), ymax=trans(lim.male)+1), fill="white") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(limits=c(0,NA), breaks=trans(yticks), labels=yticks) +
  labs(y="# genes")  + ggtitle (paste("#log2(FC)<=-", thresh.norm.FC,  " per contig - male enrichment\nContig (nb tot genes)", sep = ""))

pdf(paste("results/coverage/barplot_FC_threshold", thresh.norm.FC, "_per_sexe.pdf", sep = ""), w = 12, h = 8)
multiplot(gf, gm, cols = 1)
dev.off()
system(paste("bash src/date.sh results/coverage/barplot_FC_threshold", thresh.norm.FC, "_per_sexe.pdf", sep = ""))

##### Histogram of pvalue for FC tests at bp level
tests <- counts.genes.bp
pdf("results/coverage/pval_output_tests_genes_log2_FC.pdf", w = 7, h = 10)
par(mfrow = c(4, 2))
indicators <- c("t-test female vs male", "t-test unilateral lower", "median test FC", "t-test FC")
tests_done <- c("pval.ttest.both.a", "pval.ttest.both.lower.a", "pval.med.a", "pval.ttest.a")
for (i in seq_along(tests_done)) {
  for(type in c("BH", "BY")) {
    hist(p.adjust(tests[, tests_done[i]], type), main = paste(indicators[i], " on female vs male\nAnscombe transf. norm. counts on genes", sep = ""),
      xlab = paste("Adjusted p-value with ", type, sep = ""), breaks = 20, xaxt = "n")
    axis(side = 1, at = seq(0, 1, by = 0.05), labels = seq(0, 1, by = 0.05))
  }
}
dev.off()
system("bash src/date.sh results/coverage/pval_output_tests_genes_log2_FC.pdf")

