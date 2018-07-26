prot_busco <- read.csv('~/Documents/stage_mbelari/results/annotation/BUSCO/run_mbelari_prot/full_table_mbelari_prot.txt', sep='\t')
groups <- as.character(unique(prot_busco[prot_busco$Status=='Duplicated',1]))
info_busco <- read.csv('~/Documents/stage_mbelari/data/dataset_busco/nematoda_6231_OrthoDB9_orthogroup_info.txt', sep='\t')
info_groups <- match(groups,info_busco[,1])


