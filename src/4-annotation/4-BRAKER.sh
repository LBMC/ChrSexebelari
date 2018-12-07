export PATH=$PATH:~/Programs/augustus-3.3.1
export GENEMARK_PATH=~/Programs/gm_et_linux_64/gmes_petap/
export PATH=$PATH:~/Programs/BRAKER2/BRAKER_v2.1.0
export TOOLDIR=$PATH:~/Programs/tooldir
export PATH=$PATH:~/Programs/augustus-3.3.1/bin:~/Programs/augustus-3.3.1/scripts
export AUGUSTUS_CONFIG_PATH=~/Programs/augustus-3.3.1/config/
export AUGUSTUS_BIN_PATH=$PATH:~/Programs/augustus-3.3.1/bin
export AUGUSTUS_SCRIPTS_PATH=$PATH:~/Programs/augustus-3.3.1/scripts
BAMTOOLS_PATH=/usr/include/bamtools

BAMinput=results/annotation/STAR/Aligned.out.sorted.bam
GENOMEinput=results/annotation/RepeatMasker/final_assembly.fasta.soft.masked
OUTdir=results/annotation/BRAKER/

~/Programs/BRAKER2/BRAKER_v2.1.0/braker.pl --cores=8 --softmasking 1 --workingdir=${OUTdir} --genome=${GENOMEinput} --species=m_belari --bam=${BAMinput}


