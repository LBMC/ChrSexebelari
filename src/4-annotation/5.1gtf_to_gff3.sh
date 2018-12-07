AUGUSTUSGFFoutput=results/annotation/BRAKER/augustus.hints.gtf
GENOME=results/hybrid_assembly/DBG2OLC/final_assembly.fasta

~/Programs/BRAKER2/BRAKER_v2.1.0/gtf2gff.pl <${AUGUSTUSGFFoutput} --printExon --gff3 --out=augustus.hints.gff3


