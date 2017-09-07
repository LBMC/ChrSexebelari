  require(data.table)  
  require(dplyr)
  require(doParallel)
  require(foreach)
  require(parallel)
  require(ggplot2)
  
  ##### Compute average coverage per contig for both female and male pools
  where = "./results/coverage/"
  cov.male <- tbl_df(fread(paste(where, "2017_08_08_MRDR6_trim_Mbelari_SOFT_mapped_qual_sort.bed", sep = ""), h = F, sep = "\t", stringsAsFactors = F))
  cov.female <- tbl_df(fread(paste(where, "2017_08_08_MRDR5_trim_Mbelari_SOFT_mapped_qual_sort.bed", sep = ""), h = F, sep = "\t", stringsAsFactors = F))
  
  for (sex in c("male", "female")){
    res <- NULL
    print(sex)
    eval(parse(text = paste("cov <- cov.", sex, sep = "")))
    eval(parse(text = paste("contig <- unique(cov$V1)", sep = "")))
    interv <- round(seq(1, length(contig), length.out = 200))
    if(tail(interv, 1) != length(contig)){interv[length(interv)] <- length(contig)}
    for (i in 1:199){
      print(i)
      interv.tmp <- interv[i:(i+1)]; interv.tmp <- interv.tmp[1]:interv.tmp[2]
      for (j in interv.tmp){
        ct <- contig[j]
        tmp <- filter(cov, V1 == ct)
        s <- as.data.frame(summarise(tmp, mean=mean(V3), sd=sd(V3))); 
        V3.eff <- tmp$V3[tmp$V3>0]; 
        res <- rbind(res, data.frame(contig = ct, sex = sex, mean.depth = s$mean, sd.depth = s$sd, mean.depth.eff = mean(V3.eff), sd.depth.eff = sd(V3.eff), stringsAsFactors = F))
      }
    }
    write.table(res, paste(where, "summary_coverage_contig_inter_", sex, ".txt", sep = ""), col.names = T, quote = F, row.names = F)
  }

  ##### histogram of mean coverage: !!!!! do function !!!!!
  cov.male <- read.table(paste(where, "summary_coverage_contig_inter_male.txt", sep = ""), sep = " " , h = T) 
  pdf(paste(where, "hist_mean_coverage_over_contig_male.pdf", sep = ""))
  ggplot(data = cov.male, aes(mean.depth)) + geom_histogram(binwidth = 5, alpha = .5) + theme_bw() + geom_vline(aes(xintercept=mean(mean.depth)), linetype="dashed")+xlab("mean depth over contigs") +
  scale_x_continuous(limits=c(0, 100)) 
  dev.off()
  pdf(paste(where, "hist_mean_coverage_eff_over_contig_sex.pdf", sep = ""))
  ggplot(data = res, aes(mean.depth.eff, fill = sex)) + geom_histogram(binwidth = 5, alpha = .5) + theme_bw() + geom_vline(data = mu.df, aes(xintercept=mu.mean.depth.eff, color = sex),inherit.aes =F, linetype="dashed")+xlab("mean eff depth over contigs")
  dev.off()
  
  ##### mean cov female vs male: !!!!! do function !!!!!
  res.female <- res[which(res$sex == "female"), ]
  res.male <- res[which(res$sex == "male"), ]
  res.merge <- merge(res.female, res.male, by = "contig")
  plot(res.merge$mean.depth.x, res.merge$mean.depth.y, xlab = "mean depth over contigs (female)", ylab = "mean depth over contigs (male)")  
  
  ##### Read coverage file per contig produced by tablet
  cov.stat.female <- read.table(paste(where, "statistics_contig_tablet_MRDR5_trim_Mbelari_SOFT_mapped_qual_sort.txt", sep = ""), h = T, sep = "\t")
  cov.stat.male <- read.table(paste(where, "statistics_contig_tablet_MRDR6_trim_Mbelari_SOFT_mapped_qual_sort.txt", sep = ""), h = T, sep = "\t")
  write.table(cov.stat.male[, c("Contig", "Length")], col.names = F, row.names = F, quote = F, paste(where, "Mbelari.genome", sep = ""))
  
  ##### ABOVE IS DRAFT #####
  
  # Summary of contig length
  print(summary(cov.stat.female$Length))
 
  # Stats on library size 
  s.female <- sum(cov.stat.female$Reads)
  s.male <- sum(cov.stat.male$Reads)
  print(summary(cov.stat.female$Reads))
  print(summary(cov.stat.male$Reads))

  cov.female <- cov.stat.female$Reads/cov.stat.female$Length
  cov.male <- cov.stat.male$Reads/cov.stat.male$Length
  summary(cov.female)
  summary(cov.male)

  cov.female.3000 <- cov.stat.female[which(cov.stat.female$Length<=3000), ]$Reads/cov.stat.female[which(cov.stat.female$Length<=3000), ]$Length
  cov.male.3000 <- cov.stat.male[which(cov.stat.male$Length<=3000), ]$Reads/cov.stat.male[which(cov.stat.male$Length<=3000), ]$Length
  summary(cov.female.3000)
  summary(cov.male.3000)

  #However, the coverage is higher for small contig in male sample.
  
  mean.between.sex <- sapply(1:dim(cov.stat.female), function(x) mean(c(cov.stat.female[x, ]$Reads, cov.stat.male[x, ]$Reads)))
  diff.between.sex <- sapply(1:dim(cov.stat.female), function(x) log2(cov.stat.female[x, ]$Reads/cov.stat.male[x, ]$Reads))
  mean.between.sex.norm <- sapply(1:dim(cov.stat.female), function(x) mean(c(cov.stat.female[x, ]$Reads/s.female, cov.stat.male[x, ]$Reads/s.male)))
  diff.between.sex.norm <- sapply(1:dim(cov.stat.female), function(x) log2((cov.stat.female[x, ]$Reads/s.female)/(cov.stat.male[x, ]$Reads/s.male)))
                                  
  hist(cov.stat.female$Length) # same for male and female
  
  hist(cov.stat.female$Reads)
  hist(cov.stat.male$Reads)
  
  par(mfrow = c(4, 2))
  plot(Reads~Length, cov.stat.female, main = "female")
  plot(Reads~Length, cov.stat.male, main = "male")
  
  plot(cov.stat.female$Length, cov.stat.female$Reads/cov.stat.female$Length, main = "female", xlab = "contig length", ylab = "nb reads per contig", xlim = c(500, 3000),ylim=c(0,50))
  plot(cov.stat.male$Length, cov.stat.male$Reads/cov.stat.male$Length, main = "male", xlab = "contig length", ylab ="nb reads per contig", xlim = c(500, 3000),ylim=c(0,50))
  
  plot(cov.stat.female$Reads, cov.stat.male$Reads, xlab = "nb reads per contig female", ylab = "nb reads per contig male")
  points(cov.stat.female$Reads,cov.stat.female$Reads,col="red",type="l")
  
  plot(cov.stat.female$Length/cov.stat.female$Reads, cov.stat.male$Length/cov.stat.male$Reads, xlab = "ratio Length/nb reads per contig female", 
       ylab = "ratio Length/nb reads per contig male")
  
  plot(cov.stat.female$Reads/s.female, cov.stat.male$Reads/s.male, xlab = "nb reads per contig female/Ntot female", ylab = "nb reads per contig male/Ntot male")
  plot(cov.stat.female$Length/(cov.stat.female$Reads/s.female), cov.stat.male$Length/(cov.stat.male$Reads/s.male), xlab = "ratio Length/(nb reads per contig female/Ntot female)", 
       ylab = "ratio Length/(nb reads per contig male/Ntot male)")
  
  plot(log2((cov.stat.female$Reads/cov.stat.male$Reads)*(s.male/s.female)), xlab = "contig ID", ylab = "log2((nb reads per contig female/nb reads per contig male)*R))\nwith R = Ntot male/Ntot female")
  plot(mean.between.sex, diff.between.sex)
  
  mean.between.sex.norm <- sapply(1:dim(cov.stat.female), function(x) mean(c(log2(cov.stat.female[x, ]$Reads/s.female), log2(cov.stat.male[x, ]$Reads/s.male))))
  diff.between.sex.norm <- sapply(1:dim(cov.stat.female), function(x) log2(cov.stat.female[x, ]$Reads/s.female)-log2(cov.stat.male[x, ]$Reads/s.male))
  x11()
  plot(mean.between.sex.norm, diff.between.sex.norm)
  abline(h = 0, col = "red")
  
  require(edgeR)
  #norm.fac <- calcNormFactors(cbind(cov.stat.female$Reads, cov.stat.male$Reads), lib.size=NULL, method="TMM", logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE)
  norm=cpm(cbind(cov.stat.female$Reads, cov.stat.male$Reads), normalized.lib.sizes=T)
  cov.stat.female$cpm.Reads <- norm[, 1]
  cov.stat.male$cpm.Reads <- norm[, 2]
  
  par(mfrow = c(2,2))
  plot(cov.stat.female$Length, cov.stat.female$cpm.Reads/cov.stat.female$Length, main = "female", xlab = "contig length", ylab = "norm nb reads per contig", xlim = c(0, 5000))
  plot(cov.stat.male$Length, cov.stat.male$cpm.Reads/cov.stat.male$Length, main = "male", xlab = "contig length", ylab ="norm nb reads per contig", xlim = c(0, 5000))
 
  plot(cov.stat.female$Reads, cov.stat.male$Reads)
  plot(cov.stat.female$cpm.Reads, cov.stat.male$cpm.Reads)
  points(cov.stat.female$cpm.Reads,cov.stat.female$cpm.Reads,col="red",type="l")
  
  
  
  ##### Open bedpe files 
  require(dplyr); require(data.table)
  MRDR5 <- fread("./results/coverage/2017_08_08_MRDR5_trim_Mbelari_SOFT_mapped_qual_sort.bedpe", sep = "\t", h = F)
  MRDR5 <- tbl_df(MRDR5)
  MRDR6 <- fread("./results/coverage/2017_08_08_MRDR6_trim_Mbelari_SOFT_mapped_qual_sort.bedpe", sep = "\t", h = F)
  MRDR6 <- tbl_df(MRDR6)
  mrdr6 <- tbl_df(rbind(as.matrix(select(MRDR6, V1,V2)), as.matrix(select(MRDR6, V4,V5))))
  mrdr5 <- tbl_df(rbind(as.matrix(select(MRDR5, V1,V2)), as.matrix(select(MRDR5, V4,V5))))
  
  mrdr6=arrange(mrdr6, V1,V2); mrdr6$V2 <- as.numeric(mrdr6$V2)
  mrdr6.test = filter(mrdr6, V1 == "MBELA00001")
  
  
  
  