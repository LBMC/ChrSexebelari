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
  
##### Data file 
gff3.file.all.genome <- "data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2.gff3"
contigs.mbela <- fread("data/ReferenceGenomes/2017_09_13_Mbelari.sizes.genome", sep = "\t", h = F, stringsAsFactors = F)

compute.norm <- F
if(compute.norm){
  ##### Plot cluster size when comparing unmapped female and male reads
  clus.size <- read.csv("results/mapping/mapped/2017_11_22_cdhit.fa.clstr.size", sep = "\t", h = T, stringsAsFactors = F)
  #> dim(clus.size)
  #[1] 9079667       2
  #summary(clus.size$nb)
  #  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  #0.000    0.000    0.000    0.713    0.000 4430.000 
  pdf("results/mapping/mapped/cluster_size.pdf")
  hist(log10(clus.size$nb), breaks = 50, main = paste("log10(size of clusters) (n=", dim(clus.size)[1], ") obtained with CD-HIT\nto compare unmapped female and male reads", sep = ""), cex.main = 0.9)
  dev.off()
  system("bash src/date.sh results/mapping/mapped/cluster_size.pdf")

  ##### Pick contig name and length
  print(summary(contigs.mbela$V2))

  ##### Pick annotation to get gene regions 
  gff.mbela <- read.csv(gff3.file.all.genome, sep = "\t", h = F, stringsAsFactors = F)
  gff.genes.mbela <- gff.mbela[which(gff.mbela$V3 == "gene"), ]

  ##### Percentage of genic region:
  tot.size <- sum(contigs.mbela$V2)
  gene.size <- abs(gff.genes.mbela$V5 - gff.genes.mbela$V4)
  perc.genic <- 100*sum(gene.size)/tot.size # 52.72679
  perc.ctg <- sapply(1:dim(contigs.mbela)[1], function(x) {ctg = contigs.mbela[x]; tmp = gff.genes.mbela[which(gff.genes.mbela$V1 == ctg), ]; s = sum(abs(tmp$V5 - tmp$V4)); 100*s/contigs.mbela$V2[x]})
  #> summary(perc.ctg)
  #  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  #0.000000 0.000000 0.000000 0.000373 0.000000 6.082403

  ##### Summary nb genes per contig
  gff.genes.mbela$V10 <- sapply(gff.genes.mbela$V9, function(x) strsplit(x, "ID=")[[1]][2])
  genes.contigs.mbela <- table(gff.mbela$V1, gff.mbela$V3)
  pdf("results/coverage/number_genes_detected_in_contigs.pdf")
  barplot(table(genes.contigs.mbela[, "gene"]), main = paste("nb genes within contig (genic% = ", perc.genic, ")", sep = ""), ylab = "Frequency", xlab = "nb genes per contig")
  dev.off()

  ##### Nb of reads per contig 
  cov.female <- fread("results/coverage/2017_10_26_coverage_summary_tablet_female_JU2817_trim_Mbelari_mapped_rmdup_rg_realign_indels.txt", sep = "\t", h = T, stringsAsFactors = F)
  cov.male <- fread("results/coverage/2017_10_26_coverage_summary_tablet_male_JU2817_trim_Mbelari_mapped_rmdup_rg_realign_indels.txt", sep = "\t", h = T,  stringsAsFactors = F)
  
  ##### Summary number of reads on contigs per pool after having removed duplicates
  cov.female <- cov.female[which(cov.female$Contig %in% contigs.mbela$V1), ]
  cov.male <- cov.male[which(cov.male$Contig %in% contigs.mbela$V1), ]
  print(summary(cov.female$Reads))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#      0     108     478    2170    2455  108520
  print(summary(cov.male$Reads))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#      0      62     284    1230    1351  101118 

  ##### Mean coverage value per sample
  print(sum(cov.male$Reads*100*2)/sum(cov.male$Length)) #24.76252
  print(sum(cov.female$Reads*100*2)/sum(cov.female$Length)) #43.67202

  ##### Normalize genes counts 
  # Generate bed file for gene region coordinates
  system("bash src/convertgff_to_bed.sh data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_genes.gff3 data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_genes.bed")
  # Extract gene region in bed depth file at each bp with intersect_bed_genes.sh
  system("bash src/intersect_gene_bed_genes.sh")
 
  ##### Normalize counts at contig level, gene level and at gene bp level
  counts.contigs <- ComputeNormalizedCount(cov.female, cov.male, "contig")
  write.table(counts.contigs, "results/coverage/FC_normalized_coverage_at_contig.txt", sep= "\t", quote =F, col.names = T, row.names = F)
  system("bash src/date.sh results/coverage/FC_normalized_coverage_at_contig.txt") 

  counts.genes.female <- fread("results/coverage/2017_12_06_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_count_genes.txt", h = F, sep = "\t", stringsAsFactors = F)
  counts.genes.male <- fread("results/coverage/2017_12_06_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_count_genes.txt", h = F, sep = "\t", stringsAsFactors = F)
  counts.genes <- ComputeNormalizedCount(counts.genes.female, counts.genes.male, "gene")
  write.table(counts.genes, "results/coverage/FC_normalized_coverage_at_gene.txt", sep= "\t", quote =F, col.names = T, row.names = F)
  system("bash src/date.sh results/coverage/FC_normalized_coverage_at_gene.txt") 

  ##### Intersect genes and coverage at bp resolution
  counts.genes.female <- tbl_df(fread("results/coverage/2017_12_06_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_sort_in_genes.bed", h = F, sep = "\t", stringsAsFactors = F))
  counts.genes.male <- tbl_df(fread("results/coverage/2017_12_06_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_sort_in_genes.bed", h = F, sep = "\t", stringsAsFactors = F))
  #in_male <- setdiff(unique(counts.genes.male$V8), unique(counts.genes.female$V8))
  #in_female <- setdiff(unique(counts.genes.female$V8), unique(counts.genes.male$V8))
  #> in_female #24211 genes
  #[1] "MBELA.g12439" "MBELA.g12585" "MBELA.g14323" "MBELA.g14395" "MBELA.g14542"
  #[6] "MBELA.g17943" "MBELA.g23060"
  #> in_male #24209 genes
  #[1] "MBELA.g11678" "MBELA.g15899" "MBELA.g17246" "MBELA.g17397" "MBELA.g20649"

  ##### Compute subset of data to do merge
  counts.genes.bp <- ComputeNormalizedCount(counts.genes.female, counts.genes.male, "bp")
  write.table(counts.genes.bp, "results/coverage/FC_normalized_coverage_at_bp_within_genes.txt", sep= "\t", quote =F, col.names = T, row.names = F)
  system("bash src/date.sh results/coverage/FC_normalized_coverage_at_bp_within_genes.txt")
}

do.test.at.bp <- F
if(do.test.at.bp ){
  ##### In case of bp level, implement test per gene on estimated FC valeurs: if not gaussian, compute i) a glm based on NB to compare counts, ii) a median like based test on the log2(FC), iii) try the the ans comb transformation and perform a test on it.
  counts.bp <- tbl_df(fread("results/coverage/2017_11_30_FC_normalized_coverage_at_bp_within_genes.txt", sep= "\t", stringsAsFactors = F))
  genes <- counts.bp %>% distinct(V8); genes <- unique(genes$V8) 

  res <- NULL
  for(g in genes) {
        tmp <- filter(counts.bp, V8 == g)[, c("counts.norm.female", "counts.norm.male", "log2.norm.FC")]
	fc <- tmp$log2.norm.FC
	f <- tmp$counts.norm.female; f.a <- 2*sqrt(f+3/8)
	m <- tmp$counts.norm.male; m.a <- 2*sqrt(m+3/8)
	fc.a <- log2(f.a/m.a)
	me <- median(fc); mean <- mean(fc)
	me.f <- median(f); mean.f <- mean(f); 
	me.m <- median(m); mean.m <- mean(m);
 
	if (length(unique(fc)) == 1) {
		med.test <- NA	
		t.test <- NA
		t.test.a <- NA
		t.test.both <- NA
		t.test.both.a <- NA
	} else{
		med.test <- prop.test(sum(fc > 0), length(fc), p = 0.5, "two.sided")$p.value
		med.test.a <- prop.test(sum(fc.a > 0), length(fc.a), p = 0.5, "two.sided")$p.value		
		t.test <- t.test(fc)$p.value
		t.test.a <- t.test(fc.a)$p.value		
		t.test.both <- t.test(f, m)$p.value	
		t.test.both.a <- t.test(f.a, m.a)$p.value
		#dat <- data.frame(count = round(c(m, f)), sexe = c(rep("male", length(m)), rep("female", length(f))))
		#m1 <- MASS::glm.nb(count~sexe, data = dat, control = glm.control(maxit=200))	
		#m1.g <- glm(count~sexe, data = dat, family = "poisson")	
		#m1.q <- glm(count~sexe, data = dat, family = "quasipoisson")	
		#m1.n <- glm(count~sexe, data = dat, family = "gaussian")	
		#m0 <- MASS::glm.nb(count~1, data = dat,control = glm.control(trace = 10,maxit=200))	
		#m0.g <- glm(count~1, data = dat, family = "poisson")	
		#m0.q <- glm( count~1, data = dat, family = "quasipoisson")	
		#m0.n <- glm(count~1, data = dat, family = "gaussian")	
		#pm.anova <- anova(m0, m1)[["Pr(Chi)"]][2]
	}
	res <- rbind(res, data.frame(gene = g, median.log2.norm.FC = me, mean.log2.norm.FC = mean, median.female = me.f, median.male = me.m, mean.female = mean.f, mean.male = mean.m, pval.med = med.test, pval.med.a = med.test.a, pval.ttest = t.test, pval.ttest.a = t.test.a, pval.ttest.both = t.test.both, pval.ttest.both.a = t.test.both.a, stringsAsFactors = F))
  }
  write.table(res, "tests_FC_normalized_coverage_at_bp_within_genes.txt", sep= "\t", quote =F, col.names = T, row.names = F)
  system("bash src/date.sh results/coverage/tests_FC_normalized_coverage_at_bp_within_genes.txt")
}

##### Need to check to rerun
counts.contigs <- read.csv("results/coverage/2017_12_06_FC_normalized_coverage_at_contig.txt", sep= "\t", h = T, stringsAsFactors =  F)
counts.genes <- read.csv("results/coverage/2017_12_06_FC_normalized_coverage_at_gene.txt", sep= "\t", h = T, stringsAsFactors =  F)
counts.genes.bp <- read.csv("results/coverage/2017_11_30_tests_FC_normalized_coverage_at_bp_within_genes.txt", sep= "\t", h = T, stringsAsFactors = T)

# Plot on different scale of log2 fc and see effects of normalization
##### Histogram of log2(FC) at contig, gene and estimated at gene bp level
pdf("results/coverage/all_log2_FC.pdf", w = 10, h =5)
par(mfrow = c(1, 3))
hist(counts.contigs$log2.raw.FC, main = "Histogram of log2(female/male)) for contig", xlab = "log2(FC)", breaks = 100)
hist(counts.contigs$log2.norm.FC, breaks = 100, add = T, col = adjustcolor("red", 0.75))
legend("topleft", title = paste("at contig level (n=", dim(counts.contigs)[1], ")", sep = ""), legend = c("log2(FC) on quantile normalized counts", "log2(FC) on raw counts"), fill = c("red", "white"), cex = 0.75, bty = "n")

hist(counts.genes$log2.raw.FC, main = "Histogram of log2(female/male)) for genes", xlab = "log2(FC)", breaks = 100)
hist(counts.genes$log2.norm.FC, breaks = 100, add = T, col = adjustcolor("red", 0.75))
legend("topleft", title = paste("at gene level (n=", dim(counts.genes)[1], ")", sep = ""), legend = c("log2(FC) on quantile normalized counts", "log2(FC) on raw counts"), fill = c("red", "white"), cex = 0.75, bty = "n")

hist(counts.genes.bp$mean.log2.norm.FC, main = "Histogram of log2(female/male)) for\ngenes from bp resolution", xlab = "log2(FC)", breaks = 100,  col = adjustcolor("red", 0.75))
dev.off()

system("bash src/date.sh results/coverage/all_log2_FC.pdf")
  

