require(data.table)
source('http://bioconductor.org/biocLite.R')
biocLite('preprocessCore')
library(preprocessCore)

source('~/Documents/stage_mbelari/src/coverage_analysis/functions_claire.R')

contigs.mbela <- read.table("~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-Mbelari_hybrid_genome_sizes.txt", sep = "\t", head = F, row.names=1)
cov.female <- fread("~/Documents/stage_mbelari/results/coverage_analysis/2018_06_29_MRDR6_vs_redundans_assembly.sorted.filtered.counts", sep = "\t", h = T, stringsAsFactors = F)
cov.male <- fread("~/Documents/stage_mbelari/results/coverage_analysis/2018_06_29_MRDR6_vs_redundans_assembly.sorted.filtered.counts", sep = "\t", h = T, stringsAsFactors = F)

colnames(cov.female) <- c('Contig','Length','Reads','Unmapped')
colnames(cov.male) <- c('Contig','Length','Reads','Unmapped')

cov.female <- cov.female[which(cov.female$Contig %in% rownames(contigs.mbela)), ]
cov.male <- cov.male[which(cov.male$Contig %in% rownames(contigs.mbela)), ]

counts.contigs <- ComputeNormalizedCount(cov.female, cov.male, "contig")

counts.contigs$missing.sexe <- ""; counts.contigs$missing.sexe[which(counts.contigs$counts.raw.female == 0)] <- "female";
counts.contigs$missing.sexe[which(counts.contigs$counts.raw.male == 0)] <- "male";
write.table(counts.contigs, "2018-06-21-FC_normalized_coverage_at_contig.txt", sep = "\t", quote = F, col.names = T, row.names = F)

counts.contigs <- read.table('~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-FC_normalized_coverage_at_contig.txt', sep='\t',head=T,row.names=1)


res <- NULL
for(g in genes) {
    tmp <- as.data.frame(filter(counts.contigs, V8 == g)[, c("V1", "counts.norm.female", "counts.norm.male", "log2.norm.FC", "counts.raw.female", "counts.raw.male", "log2.raw.FC")])
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

bed.female <- fread("~/Documents/stage_mbelari/results/coverage_analysis/2018_06_21_MRDR5_vs_Hybrid_assembly.sorted.filtered.bed", sep = "\t", h = F, stringsAsFactors = F)
bed.male <- fread("~/Documents/stage_mbelari/results/coverage_analysis/2018_06_21_MRDR5_vs_Hybrid_assembly.sorted.filtered.bed", sep = "\t", h = F, stringsAsFactors = F)

contigs <- rownames(contigs.mbela)
for (c in contigs){
tmp <- bed.female[bed.female[,1]==c,]
}


##################################

library(gplots)

blast.results <- read.table("~/Documents/stage_mbelari/results/13genes/version2/2018-06-28-Genes_absent_females.txt", sep = "\t", head = F)
counts.contigs <- read.table('~/Documents/stage_mbelari/results/coverage_analysis/2018-06-21-FC_normalized_coverage_at_contig.txt', sep='\t',head=T,row.names=1)

genes.contigs <- table(blast.results[,1],blast.results[,2])
#reduce the heat map
genes.contigs <- genes.contigs[c(1,2,4,6,7,14,15),]
genes.contigs <- genes.contigs[,apply(genes.contigs,2,sum)>0]


log2fc <- counts.contigs[colnames(genes.contigs),'log2.norm.FC']
#genes.contigs <- rbind(genes.contigs,log2fc)

Contigs.low.fc <- list()

for (n in 1:nrow(genes.contigs)){
contigs.present <- which(genes.contigs[c(n),]>0)
Contigs.low.fc[[n]] <- genes.contigs['log2fc',names(contigs.present)]	
}

names(Contigs.low.fc) <- rownames(genes.contigs)


log2fc.contigs <- matrix(log2fc,nrow=nrow(genes.contigs),ncol=ncol(genes.contigs), byrow=TRUE)
log2fc.contigs[which(genes.contigs==0)]<-NA
my_col <- colorRampPalette(c("green","black","red"))(256)
rownames(log2fc.contigs) <- rownames(genes.contigs)
colnames(log2fc.contigs) <- colnames(genes.contigs)
rownames(log2fc.contigs) <- sub('MBELA\\.','',rownames(log2fc.contigs))

heatmap.2(t(log2fc.contigs), col= my_col, trace='none', Rowv=FALSE, Colv=FALSE,keysize=1,key.title=NA, key.xlab='Log2FC', cexRow=1, cexCol=1.2, margins=c(5.7,7) , notecol='white', density.info='none')
















