library(gridExtra)
require(modeest)
require(Sushi)
require(qvalue)
require(data.table)
require(foreach)
require(doParallel)
require(VGAM)
source("src/func/functions.R")

do.fisher <- F
alpha.threshold <- 0.01

if (do.fisher) {
  data.filt.SNP <- read.csv("results/call_var/2017_12_08_merge.sexe.all.filt.SNP.txt", sep = "\t", h = T, stringsAsFactors = F)

  # Fisher tests are not computed on all SNPs!!
  # Pick same alt allele for female and male
  data.filt.SNP.same.alt <- data.filt.SNP[which(data.filt.SNP$same.ALT == 1), ]

  # Compute freq male and freq female
  data.filt.SNP.same.alt$freq.alt.female <- data.filt.SNP.same.alt$count.alt.female/(data.filt.SNP.same.alt$count.ref.female+data.filt.SNP.same.alt$count.alt.female)
  data.filt.SNP.same.alt$freq.alt.male <- data.filt.SNP.same.alt$count.alt.male/(data.filt.SNP.same.alt$count.ref.male+data.filt.SNP.same.alt$count.alt.male)

  ##################
  ##### Fisher #####
  ##################

  # Add raw pvalue associated to Fisher test done on raw counts
  pval.fisher <- ComputePvalueFisher(input.df = data.filt.SNP.same.alt[, c("count.ref.male", "count.alt.male", "count.ref.female", "count.alt.female")])
  data.filt.SNP.same.alt$pval.fisher.raw.count <- pval.fisher

  # Adjust Fisher pvalue with BH
  data.filt.SNP.same.alt$pval.adjBH.fisher.raw.count <- p.adjust(data.filt.SNP.same.alt$pval.fisher.raw.count, method = "BH")

  write.table(data.filt.SNP.same.alt, "results/fisher/merge.sexe.all.filt.SNP.same.alt.fisher.txt", sep = "\t",  quote = F, row.names = F)
  system("bash src/date.sh results/fisher/merge.sexe.all.filt.SNP.same.alt.fisher.txt")

  write.table(data.filt.SNP.same.alt[which(data.filt.SNP.same.alt$genes != ""), ], "results/fisher/merge.sexe.all.filt.SNP.same.alt.fisher.genes.txt", sep = "\t",  quote = F, row.names = F)
  system("bash src/date.sh results/fisher/merge.sexe.all.filt.SNP.same.alt.fisher.genes.txt")
}

do.H0 <- F
if(do.H0){
	##### Compute H0 distribution of count tables: UPDATE: instead using Bastide et al Methods to choose the best alpha value, I use instead a vglm fit on betabinomial to get alpha and beta parameters. Version 1.02 of VGAM was used on bioinfo space computer
	data.tests <- read.csv("results/fisher/2017_12_22_merge.sexe.all.filt.SNP.same.alt.fisher.txt", sep = "\t", h = T, stringsAsFactors = F)
	data.tests.genes <- read.csv("results/fisher/2017_12_22_merge.sexe.all.filt.SNP.same.alt.fisher.genes.txt", sep = "\t", h = T, stringsAsFactors = F)

        # for all SNPs
	fit <- vglm(cbind(c(data.tests$count.ref.female, data.tests$count.alt.female), c(data.tests$count.ref.male, data.tests$count.alt.male)) ~ 1, betabinomial, trace = T)
	saveRDS(fit, "results/fisher/params_H0_beta.RData")
	system("bash src/date.sh results/fisher/params_H0_beta.RData")

        # for SNPs within genes
	fit <- vglm(cbind(c(data.tests.genes$count.ref.female, data.tests.genes$count.alt.female), c(data.tests.genes$count.ref.male, data.tests.genes$count.alt.male)) ~ 1, betabinomial, trace = T)
	saveRDS(fit, "results/fisher/params_H0_beta.genes.RData")
	system("bash src/date.sh results/fisher/params_H0_beta.genes.RData")
}

##### Compute 10 datasets under H0 with alpha obtained above
do.simu.10x <- F
if(do.simu.10x) {
	for (type in c(".genes", "")) {
		data.tests <- read.csv(paste("2017_12_22_merge.sexe.all.filt.SNP.same.alt.fisher", type, ".txt", sep = ""), sep = "\t", h = T, stringsAsFactors = F)
		readRDS(paste("2017_12_22_params_H0_beta", type, ".RData", sep = ""))
		coeff <- Coef(fit)
		alpha <- (1-coeff["rho"])/(coeff["rho"]*(1+(1-coeff["mu"])/coeff["mu"]))
		beta <- alpha*(1-coeff["mu"])/coeff["mu"]
		ind.ord <- order(data.tests$pval.fisher.raw.count)
		data.tests.simu <- data.tests[ind.ord, ]
		pH0.simu <- NULL
		for (i in 1:50){
			print(i)
			cores = detectCores()
			cl <- makeCluster(cores[1]-1) #not to overload your computer
			registerDoParallel(cl)
			pvalH0_out <- foreach(i=1:dim(data.tests.simu)[1], .combine=rbind) %dopar% {
		          	tmp <- data.tests.simu[i, ]
		
				nref <- tmp$count.ref.female + tmp$count.ref.male
		                pmref <- rbeta(1, shape1 = alpha, shape2 = beta)
				nrefmale <- rbinom(n = 1, size = nref, prob = pmref)
				nreffemale <- nref-nrefmale
		
				nalt <- tmp$count.alt.female + tmp$count.alt.male
        	        	pmalt <- rbeta(1, shape1 = alpha, shape2 = beta)
				naltmale <- rbinom(n = 1, size = nalt, prob = pmalt)
				naltfemale <- nalt-naltmale
	
 		  		data.frame(ID = tmp$ID, count.ref.male = nrefmale, count.alt.male = naltmale, count.ref.female = nreffemale, count.alt.female = naltfemale, stringsAsFactors = F)
			}
		  	stopCluster(cl)
	
			pvalH0_out <- pvalH0_out[sapply(data.tests.simu$ID, function(x) which(pvalH0_out$ID == x)), ]
		        pval.fisher.H0 <- ComputePvalueFisher(input.df = pvalH0_out[, c("count.ref.male", "count.alt.male", "count.ref.female", "count.alt.female")])
			pH0.simu <- cbind(pH0.simu, pval.fisher.H0)
		}
		saveRDS(pH0.simu, paste("simulated_count_H0", type, ".RData", sep = ""))	
		#system(paste("bash src/date.sh results/fisher/simulated_count_H0", type, ".RData", sep = ""))	
		#pH0.simu <- cbind(pH0.simu, data.frame(pval.obs.ord = data.tests.simu$pval.fisher.raw.count))
		# For the ith rank SNP, # nb simu SNP with pvalue <= pvalue obs for Fisher and compute fdr with rank
		#pemp <- apply(pH0.simu, 1, function(x) {y <- unlist(x); length(which(y[1:10] <= y[11]))})/10
		pemp <- sapply(data.tests.simu$pval.fisher.raw.count, function(x) length(which(as.vector(pH0.simu) <= x)))/10)
		fdr <- pemp/1:length(pemp)
		data.tests.simu$FDR <- fdr
		write.table(data.tests.simu, paste("merge.sexe.all.filt.SNP.same.alt.fisher.FDR", type, ".txt", sep = ""), sep = "\t", col.names = T, row.names = F)
		#system(paste("bash src/date.sh results/call_var/merge.sexe.common.all.filt.SNP.same.alt.fisher.FDR", type, ".txt", sep = ""))
	}
}	

##### Zoom on tests output
size.genome.mbelari <- fread("data/ReferenceGenomes/2017_09_13_Mbelari.sizes.genome", sep = "\t", h = F, stringsAsFactors = F)
# annotation
gff.mbelari <- read.csv("data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2.gff3", h = F, sep = "\t", stringsAsFactors = F)
# genes annotation
gff.mbelari.genes <- gff.mbelari[which(gff.mbelari$V3 == "gene"), ]
gff.mbelari.genes$length <- abs(gff.mbelari.genes$V5-gff.mbelari.genes$V4)
# output of fisher tests
tests <- read.csv("results/fisher/2017_12_08_merge.sexe.all.filt.SNP.same.alt.fisher.FDR.txt", sep = "\t", h = T, stringsAsFactors = F)

##### Nb of significative SNP results at gene and contig level
tests.within.genes <- tests[which(tests$genes != ""), ]

tab.ctg <- table(tests$FDR<alpha.threshold, tests$CHROM); tab.signif.ctg <- tab.ctg[, which(tab.ctg[2, ]>0)]; tab.signif.ctg <- tab.signif.ctg[,rev(order(tab.signif.ctg[2,]))]
tab.ctg.within.genes <- table(tests.within.genes$FDR<alpha.threshold, tests.within.genes$CHROM); tab.signif.ctg.within.genes <- tab.ctg.within.genes[, which(tab.ctg.within.genes[2, ]>0)]; tab.signif.ctg.within.genes <- tab.signif.ctg.within.genes[,rev(order(tab.signif.ctg.within.genes[2,]))]
tab.gene <- table(tests.within.genes$FDR<alpha.threshold, tests.within.genes$gene); tab.signif.gene <- tab.gene[, which(tab.gene[2, ]>0)]; tab.signif.gene <- tab.signif.gene[,rev(order(tab.signif.gene[2,]))]
pdf("results/fisher/signif_Fisher_contig_gene_level.pdf", h = 9, w = 12)
par(mfrow = c(1,3))
barplot(tab.signif.ctg[2, ], xlab = "At least 1 SNP signif. on this contig (FDR for fisher)", ylab = "#signif per contig", xaxt = "n",
main = paste("Signif. fisher tests (", sum(tab.signif.ctg[2, ]), ") per\nSNP at ", 100*alpha.threshold, "% per contig (", dim(tab.signif.ctg)[2], ")", sep = ""))
barplot(tab.signif.ctg.within.genes[2, ], xlab = "At least 1 SNP signif. on this contig within gene\n(FDR for fisher)", ylab = "#signif per contig within gene", xaxt = "n",
main = paste("Signif. fisher tests (", sum(tab.signif.ctg.within.genes[2, ]), ") per\nSNP at ", 100*alpha.threshold, "% per contig within gene (", dim(tab.signif.ctg.within.genes)[2], ")", sep = ""))
barplot(tab.signif.gene[2, ], xlab = "At least 1 SNP signif. on this gene (FDR for fisher)", ylab = "#signif per gene", xaxt = "n",
main = paste("Signif. fisher tests (", sum(tab.signif.gene[2, ]), ") per\nSNP at ", 100*alpha.threshold, "% per gene (", dim(tab.signif.gene)[2], ")", sep = ""))
dev.off()
system("bash src/date.sh results/fisher/signif_Fisher_contig_gene_level.pdf")

##### Plot FDR representation by contig with highest nb of significant SNPs per gene
rank.signif.gene <-  colnames(tab.signif.gene)
plot <- rank.signif.gene[1:10]
l <- sapply(plot, function(x) gff.mbelari.genes[which(gff.mbelari.genes$V9 == x), ]$length)
fisher.tmp <- tests.within.genes[unlist(sapply(plot, function(x) which(tests.within.genes$genes == x))), c("CHROM", "POS", "REF", "freq.alt.female", "freq.alt.male", "pval.adjBH.fisher.raw.count", "genes", "FDR")]
offset <- 0
for (g in plot) {
	fisher.tmp.g <- fisher.tmp[which(fisher.tmp$genes == g), ]
	if (g == plot[1]){
		plot(fisher.tmp.g$POS, fisher.tmp.g$FDR, col = c("blue", "red")[as.numeric(fisher.tmp.g$freq.alt.female<fisher.tmp.g$freq.alt.male)+1], xlim = c(0, sum(l)))
		offset <- offset + l[g]
		abline(v = offset, col = "gray")
	}else{
		points(fisher.tmp.g$POS+offset, fisher.tmp.g$FDR, col = c("blue", "red")[as.numeric(fisher.tmp.g$freq.alt.female<fisher.tmp.g$freq.alt.male)+1])
		offset <- offset + l[g]
		abline(v = offset, col = "gray")
	}
}

n <- 2
ctg <- size.genome.mbelari$V1[rev(order(size.genome.mbelari$V2))[1:n]]
length.ctg <- size.genome.mbelari$V2[rev(order(size.genome.mbelari$V2))[1:n]]
names(length.ctg) <- ctg
offset <- 0
for (c in ctg) {
	tmp <- tests[which(tests$CHROM == c), ]
	if (c == ctg[1]){
		plot(tmp$POS, tmp$FDR*c(-1, 1)[as.numeric(tmp$freq.alt.female<tmp$freq.alt.male)+1], col = c("blue", "red")[as.numeric(tmp$freq.alt.female<tmp$freq.alt.male)+1], ylim = c(-1, 1),xlim = c(0, sum(length.ctg)), ylab = "FDR", xlab = "pos")
		offset <- offset + length.ctg[c]
		abline(v = offset, col = "gray")
	}else{
		points(tmp$POS+offset, tmp$FDR*c(-1, 1)[as.numeric(tmp$freq.alt.female<tmp$freq.alt.male)+1], col = c("blue", "red")[as.numeric(tmp$freq.alt.female<tmp$freq.alt.male)+1])
		offset <- offset + length.ctg[c]
		abline(v = offset, col = "gray")
	}
}














