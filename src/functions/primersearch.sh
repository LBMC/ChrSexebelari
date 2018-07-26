#Query=~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta
Query=~/Documents/stage_mbelari/data/reference_genomes/ADBT01.1.fsa_nt
primers=~/Documents/stage_mbelari/results/13genes/version2/primers_males.txt
outname=males_ecoli.txt
mismatchper=20

primersearch -seqall ${Query} -infile ${primers} -mismatchpercent ${mismatchper} -outfile ${outname}
