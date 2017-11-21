  require(data.table)  
  require(dplyr)
  require(ggplot2)
  library(devtools)
  library(Biobase)
  library(preprocessCore)
  source("src/func/functions.R")
  
  ##### Pick contig name and length
  contigs.mbela <- fread("data/ReferenceGenomes/2017_09_13_Mbelari.sizes.genome", sep = "\t", h = F, stringsAsFactors = F)
  print(summary(contigs.mbela$V2))
  
  ##### Pick annotation to get gene regions 
  gff.mbela <- read.csv("data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2.gff3", sep = "\t", h = F, stringsAsFactors = F)
  gff.genes.mbela <- gff.mbela[which(gff.mbela$V3 == "gene"), ]

  ##### Summary nb genes per contig
  write.table(gff.genes.mbela, "data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_genes.gff3", sep = "\t", col.names = F, row.names = F, quote = F)
  gff.genes.mbela$V10 <- sapply(gff.genes.mbela$V9, function(x) strsplit(x, "ID=")[[1]][2])
  genes.contigs.mbela <- table(gff.mbela$V1, gff.mbela$V3)
  pdf("results/coverage/number_genes_detected_in_contigs.pdf")
  barplot(table(genes.contigs.mbela[, "gene"]), ylab = "Frequency", xlab = "nb genes per contig")
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

  ##### Normalize genes counts at bp level
  # Generate bed file for gene region coordinates
  system("bash src/convertgff_to_bed.sh data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_genes.gff3 data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_genes.bed")
  # Extract gene region in bed depth file at each bp with intersect_bed_genes.sh
  system("bash src/intersect_gene_bed_genes.sh")
 
  ##### Normalize genes counts at bp level within contig, within gene and at gene bp level
  counts.contigs <- ComputeNormalizedCount(cov.female, cov.male, "contig")

  counts.genes.female <- fread("results/coverage/2017_11_12_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_count_genes.txt", h = F, sep = "\t", stringsAsFactors = F)
  counts.genes.male <- fread("results/coverage/2017_11_12_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_count_genes.txt", h = F, sep = "\t", stringsAsFactors = F)
  counts.genes <- ComputeNormalizedCount(counts.genes.female, counts.genes.male, "gene")

  counts.genes.female <- fread("results/coverage/2017_11_18_MRDR5_trim_Mbelari_mapped_rmdup_rg_realign_indels_sort_in_genes.bed", h = F, sep = "\t", stringsAsFactors = F)
  counts.genes.male <- fread("results/coverage/2017_11_18_MRDR6_trim_Mbelari_mapped_rmdup_rg_realign_indels_sort_in_genes.bed", h = F, sep = "\t", stringsAsFactors = F)
  counts.genes.bp <- ComputeNormalizedCount(counts.genes.female, counts.genes.male, "bp")
  write.table(counts.genes.bp, "results/coverage/2017_11_18_FC_normalized_coverage_at_bp_within_genes.txt", sep= "\t", quote =F, col.names = T, row.names = F)

  ##### Histogram of log2(FC)
  pdf("results/coverage/2017_11_12_all_log2_FC.pdf")
  hist(counts.genes$log2.raw.FC, main = "Histogram of log2(female/male))", xlab = "log2(FC)", breaks = 100)
  hist(counts.genes$log2.norm.FC, main = "Histogram of log2(norm female/norm male)", xlab = "log2(norm FC)", breaks = 100, add = T, col = adjustcolor("gray", 0.75))
  legend("toleft", c("log2(FC) on quantile normalized counts", "log2(FC) on raw counts"), fill = c("gray", "white"), cex = 0.75, bty = "n")
  dev.off()

  ##### Zoom on genes with abs(log2(FC))>1
  signift <- counts.genes[which(abs(counts.genes$log2.norm.FC)>1), ]
  max.log2.norm.FC.per.ctg <- unique(signift$V1)[rev(order(sapply(unique(signift$V1), function(x) max(abs(signift[which(signift$V1 == x), ]$log2.norm.FC)))))]
  tab <- table(signift$V1)
  coord.max.log2.norm.FC.per.ctg <- matrix(seq_along(max.log2.norm.FC.per.ctg ), nrow = 1, dimnames = list(NULL, max.log2.norm.FC.per.ctg ))
  n <- length(coord.max.log2.norm.FC.per.ctg)
  col <- sample(rainbow(n))
  col.max.log2.norm.FC.per.ctg <- matrix(col, nrow = 1, dimnames = list(NULL, max.log2.norm.FC.per.ctg))
  pdf("results/coverage/abs_log2_norm_FC_higher_than_1.pdf", w = 16)
  plot(coord.max.log2.norm.FC.per.ctg[1, as.character(signift$V1)], signift$log2.norm.FC, pch = as.numeric(signift$V1)%%4, col = col.max.log2.norm.FC.per.ctg[1, as.character(signift$V1)], xlab = "Ordered contig by max(gfold value per contig) for abs(log2.norm.FC)>1\n1 color = 1 contig", ylab = "log2.norm.FC", cex = 0.75)
  abline(h=0, col = "gray")
  text(20,1, paste( length(which(signift$log2.norm.FC>1)), " genes on ", length(unique(signift$V1[which(signift$log2.norm.FC>1)])), " contigs", sep = ""))
  text(20,-1, paste( length(which(signift$log2.norm.FC< -1)), " genes on ", length(unique(signift$V1[which(signift$log2.norm.FC < -1)])), " contigs", sep = ""))
  dev.off()

  

