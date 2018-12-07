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

##### Quantile normalisation #####
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
   
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
   
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

##### ComputeNormalizedCount #####
ComputeNormalizedCount <- function(counts.female, counts.male, level = "gene"){
  if(level == "bp") {
    #counts.genes <- data.table:::merge.data.table(counts.female[, paste("V", c(8, 10,1, 6,7,2,4), sep = "")], counts.male[, paste("V", c(8,2,4), sep = "")], by = c("V8", "V4"))
    #counts.genes <- join(counts.female[, paste("V", c(1, 8, 10,1, 6,7,2,4), sep = "")], counts.male[, paste("V", c(8,2,4), sep = "")], by = c("V8", "V4"))
    #colnames(counts.genes)[7] <- "counts.raw.female"; colnames(counts.genes)[8] <- "counts.raw.male"

    distinct_gene_male = counts.male %>% distinct(V8); distinct_gene_male <- unique(distinct_gene_male$V8) 
    distinct_gene_female = counts.female %>% distinct(V8); distinct_gene_female <- unique(distinct_gene_female$V8) 
    all_gene <- unique(intersect(distinct_gene_male, distinct_gene_female) )
    
    # genes in common
    #tmp_male <- select(filter(counts.male, V8 %in% all_gene), "V4")
    #tmp_female <- select(filter(counts.female, V8 %in% all_gene), paste("V", c(1, 8, 10, 6,7,2,4), sep = ""))	
    tmp_male <- filter(counts.male, V8 %in% all_gene)[, "V4"]
    tmp_female <- filter(counts.female, V8 %in% all_gene)[, paste("V", c(1, 8, 10, 6,7,2,4), sep = "")]
    counts.g <- cbind(tmp_female, tmp_male)	
    colnames(counts.g)[7] <- "counts.raw.female"; colnames(counts.g)[8] <- "counts.raw.male"

    # genes in female
    in_female <- setdiff(unique(counts.female$V8), unique(counts.male$V8))
    tmp_female <- filter(counts.female, V8 %in% in_female)[, paste("V", c(1, 8, 10, 6,7,2,4), sep = "")]
    colnames(tmp_female)[7] <- "counts.raw.female"; 
    counts.g <- rbind(counts.g, cbind(tmp_female, data.frame(counts.raw.male = rep(0, dim(tmp_female)[1]), stringsAsFactors = F)))

    # genes in male
    in_male <- setdiff(unique(counts.male$V8), unique(counts.female$V8))
    tmp_male <- filter(counts.male, V8 %in% in_male)[, paste("V", c(1, 8, 10, 6,7,2,4), sep = "")]
    keep <- tmp_male[,7]
    tmp_male[, 7] <- 0	
    colnames(tmp_male)[7] <- "counts.raw.female"; 
    counts.g <- rbind(counts.g, cbind(tmp_male, data.frame(counts.raw.male = keep$V4, stringsAsFactors = F)))
    counts.genes <- counts.g
    rm(counts.g)    
    rm(counts.male)    
    rm(counts.female)
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

# Multiple plot function: put several ggplot graphs saved in an object in the same 
# layout: from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

