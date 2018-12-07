AUGUSTUSGFFoutput=results/annotation/BRAKER/augustus.hints.gff
GENOME=results/hybrid_assembly/DBG2OLC/final_assembly.fasta

~/Programs/BRAKER2/BRAKER_v2.1.0/getAnnoFasta.pl ${AUGUSTUSGFFoutput} --seqfile=${GENOME}

