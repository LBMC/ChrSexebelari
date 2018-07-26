GENOMEinput=results/annotation/RepeatMasker/final_assembly.fasta.soft.masked
GENOMEoutdir=results/annotation/STAR
N=
RNASeq1=
RNASeq2=

STAR --runThreadN ${N} --runMode genomeGenerate --genomeDir ${GENOMEoutdir} --genomeFastaFiles ${GENOMEinput}

STAR --runThreadN ${N} --genomeDir ${GENOMEoutdir} --readFilesIn ${RNASeq1} ${RNASeq2}
