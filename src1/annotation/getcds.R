library(seqinr)
augustus <- read.table('augustus.hints.gff', head=F, sep='\t')
transcripts_id  <- as.character(augustus[augustus[,3]=='transcript',9])
gene_id  <- as.character(augustus[augustus[,3]=='gene',9])

augustus_seqs <- readLines('cds_prot.txt')
cds <- augustus_seqs[grep('[atgc]{5,}',augustus_seqs)]
starts <- grep('coding sequence',cds)
starts <- c(starts,(length(cds)+1))

CDS_seqs <- list()

for (n in 1:length(transcripts_id)){
	tmp <- cds[starts[n]:(starts[n+1]-1)]
	tmp <- sub('\\# coding sequence \\= \\[','',tmp)
	tmp <- sub('\\# ','', tmp)
	tmp <- sub('\\]','', tmp)
	tmp2c <- c()

	for (x in 1:length(tmp)){
		tmp2c <- c(tmp2c,s2c(tmp[x]))
	}
	CDS_seqs[[n]] <- tmp2c

}

write.fasta(CDS_seqs,transcripts_id, 'mbelari_hybrid_cds.fasta')


#27588 genes
#30020 transcripts


