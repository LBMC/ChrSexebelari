/usr/local/bin/FastQC/fastqc --outdir results/mapping/mapped/fastqc results/mapping/mapped/2017_08_08_MRDR6_trim_Mbelari_SOFT_end2end_mapped_sort.fastq

/usr/local/bin/FastQC/fastqc --outdir results/mapping/mapped/fastqc results/mapping/mapped/2017_08_08_MRDR5_trim_Mbelari_SOFT_end2end_mapped_sort.fastq

/usr/local/bin/multiqc results/mapping/mapped/fastqc/*_fastqc.zip -o results/mapping/mapped/multiqc

