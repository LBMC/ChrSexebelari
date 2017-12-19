require(data.table)

# Convert 0-based bed obtained after gff3 to bed transformation in 1-based
# systeme coordinates
bed <- fread("data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_genes.bed", sep = "\t", 
  h = F, stringsAsFactors = F)
bed$V2 <- bed$V2 + 1
write.table(bed, "data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_genes_1based.bed", 
  sep = "\t", col.names = F, row.names = F, quote = F)
system("bash src/date.sh data/ReferenceGenomes/Mesorhabditis_belari_JU2817_v2_genes_1based.bed")
