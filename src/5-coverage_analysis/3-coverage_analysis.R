require(data.table)
require(dplyr)
library(preprocessCore)
source("~/Documents/stage_mbelari/src/coverage_analysis/functions_claire.R")

#Extract gff3 at gene level

gff.mbela <- read.csv("~/Documents/stage_mbelari/results/annotation/BRAKER/augustus.hints.gff3", 
  sep = "\t", h = F, stringsAsFactors = F)
gff.genes.mbela <- gff.mbela[which(gff.mbela$V3 == "gene"), ]
write.table(gff.genes.mbela, "~/Documents/stage_mbelari/results/coverage_analysis/Mesorhabditis_belari_JU2817_hybrid_assembly_genes.gff3", 
  col.names = F, row.names = F, quote = F, sep = "\t")

contigs.mbela <- fread("~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt", sep = "\t", 
  h = F, stringsAsFactors = F)

#the gff3 file is converted to bed
system("bash ~/Documents/stage_mbelari/src/5-coverage_analysis/convertgff_to_bed.sh ~/Documents/stage_mbelari/results/coverage_analysis/Mesorhabditis_belari_JU2817_hybrid_assembly_genes.gff3 ~/Documents/stage_mbelari/results/coverage_analysis/Mesorhabditis_belari_JU2817_hybrid_assembly_genes.bed")

#the previous output bed is formated 
source("~/Documents/stage_mbelari/src/5-coverage_analysis/format_bed.R")

#
system("bash ~/Documents/stage_mbelari/src/5-coverage_analysis/intersect_bed_genes.sh")
system("bash ~/Documents/stage_mbelari/src/5-coverage_analysis/get_counts_per_gene.sh")

#Contig level

#outputs from 1-bam_read_count.sh
cov.female <- fread("~/Documents/stage_mbelari/results/coverage_analysis/2018-06-12-MRDR5_vs_Hybrid_assembly.sorted.filtered.counts", sep = "\t", h = T, stringsAsFactors = F)
cov.male <- fread("~/Documents/stage_mbelari/results/coverage_analysis/2018-06-12-MRDR6_vs_Hybrid_assembly.sorted.filtered.counts", sep = "\t", h = T, stringsAsFactors = F)

colnames(cov.female) <- c('Contig','Length','Reads','Unmapped')
colnames(cov.male) <- c('Contig','Length','Reads','Unmapped')

cov.female <- cov.female[which(cov.female$Contig %in% rownames(contigs.mbela)), ]
cov.male <- cov.male[which(cov.male$Contig %in% rownames(contigs.mbela)), ]

counts.contigs <- ComputeNormalizedCount(cov.female, cov.male, "contig")

counts.contigs$missing.sexe <- ""; counts.contigs$missing.sexe[which(counts.contigs$counts.raw.female == 0)] <- "female";
counts.contigs$missing.sexe[which(counts.contigs$counts.raw.male == 0)] <- "male";
write.table(counts.contigs, "~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-FC_normalized_coverage_at_contig.txt", sep = "\t", quote = F, col.names = T, row.names = F)

counts.contigs <- read.table('~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-FC_normalized_coverage_at_contig.txt', sep='\t',head=T,row.names=1)


#Gene level
counts.genes.female <- fread("~/Documents/stage_mbelari/results/coverage_analysis/2018-07-24-MRDR5_vs_Hybrid_assembly.sorted.count.genes.txt", 
    h = F, sep = "\t", stringsAsFactors = F)
counts.genes.male <- fread("~/Documents/stage_mbelari/results/coverage_analysis/2018-07-24-MRDR6_vs_Hybrid_assembly.sorted.count.genes.txt", 
    h = F, sep = "\t", stringsAsFactors = F)
counts.genes <- ComputeNormalizedCount(counts.genes.female, counts.genes.male, "gene")
counts.genes$missing.sexe <- ""; counts.genes$missing.sexe[which(counts.genes$counts.raw.female == 0)] <- "female"; 
counts.genes$missing.sexe[which(counts.genes$counts.raw.male == 0)] <- "male"; 
write.table(counts.genes, "~/Documents/stage_mbelari/results/coverage_analysis/2018-07-24-FC_normalized_coverage_at_gene.txt", 
    sep = "\t", quote = F, col.names = T, row.names = F)

#At bp level
 counts.genes.female <- tbl_df(fread("~/Documents/stage_mbelari/results/coverage_analysis/2018-07-17-MRDR5_vs_Hybrid_assembly_sort_in_genes_bp.bed", 
    h = F, sep = "\t", stringsAsFactors = F))
  counts.genes.male <- tbl_df(fread("~/Documents/stage_mbelari/results/coverage_analysis/2018-07-17-MRDR6_vs_Hybrid_assembly_sort_in_genes_bp.bed", 
    h = F, sep = "\t", stringsAsFactors = F))

 counts.genes.bp <- ComputeNormalizedCount(counts.genes.female, counts.genes.male, 
    "bp")
  write.table(counts.genes, "~/Documents/stage_mbelari/results/coverage_analysis/2018-07-25-FC_normalized_coverage_at_bp_within_genes.txt", 
    sep = "\t", quote = F, col.names = T, row.names = F)

#Analysis
 gene <-read.table('~/Documents/stage_mbelari/results/coverage_analysis/2018-07-24-FC_normalized_coverage_at_gene.txt',
sep="\t", header=T)

#genes and contigs without counts

counts.contigs <- read.csv("~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-FC_normalized_coverage_at_contig.txt",  sep = "\t", h = T, stringsAsFactors = F)
counts.genes <- read.csv("~/Documents/stage_mbelari/results/coverage_analysis/2018-07-24-FC_normalized_coverage_at_gene.txt",  sep = "\t", h = T, stringsAsFactors = F)

counts.genes.missing.one.sexe <- counts.genes[which(counts.genes$missing.sexe != ""), ]
counts.genes.missing.one.sexe <- counts.genes.missing.one.sexe[order(counts.genes.missing.one.sexe$missing.sexe), ]
counts.contigs.missing.one.sexe <- counts.contigs[which(counts.contigs$missing.sexe != ""), ]
counts.contigs.missing.one.sexe <- counts.contigs.missing.one.sexe[order(counts.contigs.missing.one.sexe$missing.sexe), ]
write.table(counts.genes.missing.one.sexe, "~/Documents/stage_mbelari/results/coverage_analysis/2018-07-24-counts_per_genes_raw_counts_not_present_in_both_sexe.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(counts.contigs.missing.one.sexe, "~/Documents/stage_mbelari/results/coverage_analysis/2018-07-24-counts_per_contigs_raw_counts_not_present_in_both_sexe.txt", sep = "\t", col.names = T, row.names = F, quote = F)

genes.females.absent <- counts.genes.missing.one.sexe[counts.genes.missing.one.sexe$missing.sexe=='female',1]
write.table(genes.females.absent,'~/Documents/stage_mbelari/results/coverage_analysis/2018-07-24-Genes_absent_in_females.txt', col.names=F,row.names=F, quote=F)



#tests

counts.bp <- tbl_df(fread("~/Documents/stage_mbelari/results/coverage_analysis/2018-07-25-FC_normalized_coverage_at_bp_within_genes.txt", 
    sep = "\t", stringsAsFactors = F))

  counts.genes.missing.one.sexe <- read.csv("~/Documents/stage_mbelari/results/coverage_analysis/2018-07-24-counts_per_genes_raw_counts_not_present_in_both_sexe.txt", sep = "\t", h = T, stringsAsFactors = F)
  counts.contigs.missing.one.sexe <- read.csv("~/Documents/stage_mbelari/results/coverage_analysis/2018-07-24-counts_per_contigs_raw_counts_not_present_in_both_sexe.txt", sep = "\t", h = T, stringsAsFactors = F)

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
  counts.genes <- read.csv("~/Documents/stage_mbelari/results/coverage_analysis/2018-07-24-FC_normalized_coverage_at_gene.txt", 
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

  write.table(res, "~/Documents/stage_mbelari/results/coverage_analysis/2018-07-25-tests_FC_normalized_coverage_at_bp_within_genes.txt", sep= "\t", quote =F, col.names = T, row.names = F)

counts.genes.bp <- read.csv("~/Documents/stage_mbelari/results/coverage_analysis/2018-07-25-tests_FC_normalized_coverage_at_bp_within_genes.txt", 
  sep = "\t", h = T, stringsAsFactors = T)
counts.genes.bp.raw <- counts.genes.bp[which(counts.genes.bp$type == "raw"), ]
counts.genes.bp <- counts.genes.bp[which(counts.genes.bp$type == "norm"), ]


##### Investigation on features missing in one sexe

pdf("2018-07-25-counts_per_features_raw_counts_not_present_in_both_sexe.pdf")
par(mfrow = c(1, 2))
barplot(table(counts.genes.missing.one.sexe$missing.sexe), main = "#genes missing in one sexe", xlab = "sexe")
barplot(table(counts.contigs.missing.one.sexe$missing.sexe), main = "#contigs missing in one sexe", xlab = "sexe")
dev.off()


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

pdf("2018-07-25-counts_per_features_not_present_in_both_sexe.pdf", w = 9, h = 4)
multiplot(graw, gnorm, craw, cnorm, cols = 2)
dev.off()

##### Histogram of log2(FC) at contig, gene and estimated at gene bp level
pdf("2018-07-25-all_log2_FC_2.pdf", w = 10, h = 5)
par(mfrow = c(1, 2))
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
#hist(counts.genes.bp.raw$mean.fc, main = "Histogram of log2(female/male)) for\ngenes from bp resolution", 
#  xlab = "log2(FC)", breaks = 100)
#hist(counts.genes.bp$mean.fc, breaks = 100, add = T, col = adjustcolor("red", 0.75))
#legend("topleft", title = paste("at gene bp level (n=", dim(counts.genes)[1], ")", sep = ""), 
#  legend = c("log2(FC) on quantile normalized counts", "log2(FC) on raw counts"), 
#  fill = c("red", "white"), cex = 0.75, bty = "n")
dev.off()

##### Zoom on data with abs(FC)>=threshold 
thresh.norm.FC <- 1
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


##### Histogram of pvalue for FC tests at bp level
tests <- counts.genes.bp
pdf("2018-07-25-pval_output_tests_genes_log2_FC.pdf", w = 7, h = 10)
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

#######
s.counts.genes.bp <- counts.genes.bp[,c(1,2,4,5,7,8,10,11,19,20,23,25,27,29,32)]
s.counts.genes.bp <- data.frame(s.counts.genes.bp,p.adjust(s.counts.genes.bp$pval.med.a, 'BH'),p.adjust(s.counts.genes.bp$pval.ttest.a, 'BH'),p.adjust(s.counts.genes.bp$pval.ttest.both.a, 'BH'),p.adjust(s.counts.genes.bp$pval.ttest.both.lower.a, 'BH'))
colnames(s.counts.genes.bp)[16:19] <- c('qval.med.a','qval.ttest.a','qval.ttest.both.a','qval.ttest.both.lower.a')

#Take zoom on female absent genes
absent.fem <- s.counts.genes.bp[s.counts.genes.bp$is.gene.missing.sexe=='female',]
absent.fem <- absent.fem[which((absent.fem$qval.med.a!='NA' & absent.fem$qval.ttest.a!='NA' & absent.fem$qval.ttest.both.a!='NA' & absent.fem$qval.ttest.both.lower.a!='NA')==TRUE),]
absent.fem[]

head(s.counts.genes.bp[order(s.counts.genes.bp$pval.ttest.both.lower.a, decreasing=F),])











