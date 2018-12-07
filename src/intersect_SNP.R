#!/usr/bin/Rscript
rm(list = ls())
library("tidyverse")
library("seqinr")

args <- c(
  "results/SNP/vcf_samtools/normal_sample_raw.csv",
  "results/SNP/vcf_samtools/tumor_sample_filtered.csv",
  "results/fasta/DBG2OLC_output2_filtered.fasta",
  "data/list_of_enzymes.csv"
  )
seq_restric_size <- 21

args <- commandArgs(trailingOnly = TRUE)
snp_a <- read_delim(args[1], delim = "\t") %>%
  mutate(cords = paste0(CHROM, POS))
snp_b <- read_delim(args[2], delim = "\t") %>%
  mutate(cords = paste0(CHROM, POS))

only_b <- snp_b %>%
  select(cords) %>%
  setdiff((snp_a %>% select(cords))) %>%
  pull(cords)

snp <- snp_b %>%
  filter(cords %in% only_b) %>%
  filter(tumor_sample.GT %in% c("A/A", "T/T", "G/G", "C/C")) %>%
  mutate(REF = do.call(rbind, strsplit(AD, split = ",", fixed = TRUE))[,1],
         VAR = do.call(rbind, strsplit(AD, split = ",", fixed = TRUE))[,2],
         REF = as.integer(REF),
         VAR = as.integer(VAR),
         tumor_sample.AD = NULL,
         cords = NULL
  ) %>%
  filter(REF == 0) %>%
  filter(VAR >= 10) %>%
  arrange(CHROM, desc(VAR))

fastafile <- read.fasta(file = args[3],
                        as.string = TRUE)

snp$seq_list <- snp %>%
  apply(1, FUN = function(x, fastafile, POS, CHROM, seq_restric_size){
      begin <- as.integer(x[ POS ]) - ((seq_restric_size - 1) / 2)
      end <- as.integer(x[ POS ]) + ((seq_restric_size - 1) / 2)
      chrom <- x[ CHROM ]
      seq_restric <- fastafile[[ chrom ]] %>%
        c2s() %>%
        substr(begin, end) %>%
        s2c()
      seq_restric[((seq_restric_size - 1) / 2) + 1] <-
        toupper(seq_restric[((seq_restric_size - 1) / 2) + 1])
      print(paste0(chrom, ":", begin, "-", end, " ", seq_restric %>% c2s()))
      return(seq_restric %>% c2s())
    },
    fastafile = fastafile,
    POS = which(colnames(snp) %in% "POS"),
    CHROM = which(colnames(snp) %in% "CHROM"),
    seq_restric_size = seq_restric_size
  )

snp %>%
  write_csv(paste0(args[2], "only.csv" ))


snp <- read_csv(paste0(args[2], "only.csv" ))
enzyme_list <- read_csv(args[4]) %>%
  mutate(size = nchar(seq))

snp <- snp %>%
  mutate(enzyme = NA,
         enzyme_pos = NA,
         enzyme_seq = NA)
for (i in seq_len(nrow(enzyme_list))) {
  enzyme <- enzyme_list[i, ]
  enzyme_search <- gregexpr(toupper(enzyme$seq),
                            toupper(snp$seq_list),
                            fixed = TRUE) %>%
    lapply(FUN = function(x, enzyme, seq_restric_size) {
      if (x[1] < 0) {
        return(c(NA, NA, NA))
      } else {
        if (x[1] <= (seq_restric_size - 1) / 2 + 1 &
            x[1] + enzyme$size >=  (seq_restric_size - 1) / 2 + 1) {
        return(c(enzyme$enzyme, x[1], enzyme$seq))
        }
        return(c(NA, NA, NA))
      }
    },
    enzyme = enzyme,
    seq_restric_size = seq_restric_size)
  enzyme_search <- do.call(rbind, enzyme_search)
  snp[is.na(snp$enzyme), c("enzyme", "enzyme_pos", "enzyme_seq")] <-
    enzyme_search[is.na(snp$enzyme), ]
}

snp %>%
  filter(!is.na(enzyme)) %>%
  write_csv(paste0(args[2], "only_enzyme.csv" ))

