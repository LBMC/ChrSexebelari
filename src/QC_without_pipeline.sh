# QC and unmapped data processing
# done on local computer 

fastqc --quiet --threads 1 --outdir results/quality_control/unmapped/ results/mapping/unmapped/2017_10_30_MRDR6*_unmapped.fastq &&\
fastqc --quiet --threads 1 --outdir results/quality_control/unmapped/ results/mapping/unmapped/2017_10_30_MRDR5*_unmapped.fastq &&\

multiqc results/quality_control/unmapped/ results/quality_control/unmapped/
