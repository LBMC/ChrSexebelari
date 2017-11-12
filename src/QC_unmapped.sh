fastqc --quiet --threads 1 --outdir results/quality_control/unmapped/ results/mapping/2017_10_30_MRDR6*_unmapped.fastq &&\
fastqc --quiet --threads 1 --outdir results/quality_control/unmapped/ results/mapping/2017_10_30_MRDR5*_unmapped.fastq &&\

multiqc results/quality_control/unmapped/ results/quality_control/unmapped/
