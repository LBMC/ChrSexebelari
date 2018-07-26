GENOMEinput=results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta
databaserepeats=results/dnaPipeTE/mbelari_illumina550/Trinity.fasta
outputdir=results/annotation/RepeatMasker

~/Programs/RepeatMasker/RepeatMasker -xsmall -lib ${databaserepeats} ${GENOMEinput} -dir ${outputdir} 
#remove -xsmall for hard masking
