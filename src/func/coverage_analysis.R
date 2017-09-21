  require(data.table)  
  require(dplyr)
  require(doParallel)
  require(foreach)
  require(parallel)
  require(ggplot2)
  
  ##### Compute average coverage per contig for both female and male pools
  contigs.mbela <- read.table("data/ReferenceGenomes/2017_09_12_contigs_short_name_Mbelari.txt", h = F, sep = "\t", stringsAsFactors = F)
  
  cov.female <- fread("results/coverage/2017_09_14_coverage_summary_tablet_female_JU2817_trim_Mbelari_mapped_sort.txt", sep = "\t", h = T, stringsAsFactors = F)
  cov.male <- fread("results/coverage/2017_09_14_coverage_summary_tablet_male_JU2817_trim_Mbelari_mapped_sort.txt", sep = "\t", h = T,  stringsAsFactors = F)
  
  cov.female <- cov.female[which(cov.female$Contig %in% contigs.mbela$V1), ]
  cov.male <- cov.male[which(cov.male$Contig %in% contigs.mbela$V1), ]
  print(summary(cov.female$Reads))
  print(summary(cov.male$Reads))

  # Mean coverage value per sample
  print(sum(cov.male$Reads*100)/sum(cov.male$Length))
  print(sum(cov.female$Reads*100)/sum(cov.female$Length))

  # Define 10 longest contigs  
  longer_10contigs <- cov.female$Contig[rev(order(cov.female$Length))[1:10]]
	
  # Open bed file and keep Mbelari contigs (Eclo contigs have 0 read consistently to previous filter)
  cov.male.bed <- tbl_df(fread("results/coverage/2017_09_13_MRDR6_trim_Mbelari_mapped_sort.bed", h = F, sep = "\t", stringsAsFactors = F))
  cov.male.bed <- filter(cov.male.bed, V1 %in% contigs.mbela$V1)

  cov.female.bed <- tbl_df(fread("results/coverage/2017_09_13_MRDR5_trim_Mbelari_mapped_sort.bed", h = F, sep = "\t", stringsAsFactors = F))
  cov.female.bed <- filter(cov.female.bed, V1 %in% longer_10contigs)

  pdf("results/coverage/raw_coverage_per_contig.pdf", w = 7, h = 4)
  par(mfrow = c(1, 2))
  hist(cov.female$Reads, main = "female", xlab = "coverage per contig")	
  hist(cov.male$Reads, main = "male", xlab = "coverage per contig")
  dev.off()

 pdf("results/coverage/raw_coverage_per_bp_contig.pdf", w = 7, h = 4)
  m <- sample(1:dim(cov.female.bed)[1], 100000)
  par(mfrow = c(1, 2))
  hist(cov.female.bed$V3[m], main = "female", xlab = "100,000 coverage per bp per contig")	
  hist(cov.male.bed$V3[m], main = "male", xlab = "100,000 coverage per per bp contig")
  dev.off()

  # Ratio of normalized coverage per contig
  cov <- merge(cov.female, cov.male, by = "Contig")
  cov$ratio <- log2((cov$Reads.x/497)/(cov$Reads.y/295))
  cov$contig.end <- cumsum(cov$Length.x)
  cov$contig.start <- c(1, cov$contig.end[1:(dim(cov)[1]-1)] + 1)
  #cov.circos <- cov[, c("Contig", "contig.start", "contig.end", "Contig", "ratio")]
  #colnames(cov.circos) <- c("Chromosome", "ChromStart", "ChromEnd", "Band", "Stain")
  #chr.exclude <- NULL;
  #cyto.info <- cov.circos;
  #tracks.inside <- 1;
  #tracks.outside <- 0;
  #RCircos.Set.Core.Components(cov.circos, chr.exclude, tracks.inside, tracks.outside);
	
 pdf("results/coverage/norm_by_med_coverage_per_contig.pdf")
  plot(cov$Length.x, cov$ratio, ylab = "log2( (nb reads female/median female) / (nb reads male/median male) )", xlab = "contig length")
  dev.off()


  test.male = filter(cov.male.bed, V1 == "MBELA.09799")
  test.female = filter(cov.female.bed, V1 == "MBELA.09799")
 pdf("results/coverage/pos_norm_coverage_contigMBELA.09799.pdf")
  plot(test.male$V2, test.female$V3/295, ylab = "(nb reads/median)", main = "test on MBELA.09799 with ratio=2.01", col = "red",type = "l")
  points(test.male$V2, test.male$V3/497, col = "blue", type = "l")
  legend("topright", c("male", "female"), lty = c(1,1), col = c(1,"red"))
  dev.off()
  
  

# Plot coverage vs length ordered according to higher to lower difference 
  
