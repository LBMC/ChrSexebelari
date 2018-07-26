GENOMEinput=~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta
database=hybrid_v1

~/Programs/RepeatModeler-open-1.0.11/BuildDatabase -name ${database} -engine ncbi ${GENOMEinput}

~/Programs/RepeatModeler-open-1.0.11/RepeatModeler -engine ncbi -pa 4 -database ${database}
