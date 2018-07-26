GENOMEinput=~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta
databaserepeats=~/Documents/stage_mbelari/results/dnaPipeTE/mbelari_illumina550/Trinity.fasta
outputdir=~/Documents/stage_mbelari/results/annotation/RepeatMasker

~/Programs/RepeatMasker/RepeatMasker -xsmall -lib ${databaserepeats} ${GENOMEinput} -dir ${outputdir} 
