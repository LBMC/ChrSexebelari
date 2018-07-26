GFFfile=~/Documents/stage_mbelari/results/annotation/BRAKER/soft_masked/augustus.hints.gff3
BAMfile=~/Documents/stage_mbelari/results/13genes/2018-06-12-MRDR5_vs_Hybrid_assembly.sorted.bam

htseq-count --format=bam -i Parent ${BAMfile} ${GFFfile}


