require(data.table)

# Convert 0-based bed obtained after gff3 to bed transformation in 1-based
# systeme coordinates
bed <- fread("results/annotation/coverage_analysis/Mesorhabditis_belari_JU2817_hybrid_assembly_genes.bed", sep = "\t",
  h = F, stringsAsFactors = F)
bed$V2 <- bed$V2 + 1
bed$V3 <- bed$V3 + 1
write.table(bed, "results/annotation/coverage_analysis/Mesorhabditis_belari_JU2817_hybrid_assembly_genes_1based.bed",
  sep = "\t", col.names = F, row.names = F, quote = F)

