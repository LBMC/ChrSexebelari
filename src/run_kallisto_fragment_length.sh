# Estimate fragment length distribution on processed reads. This was done on uncompressed ref fasta file.
# done on local computer 

bin/kallisto index -i data/ReferenceGenomes/m_belari_kall_index data/ReferenceGenomes/mesorhabditis_belari_assembly.2.0.SOFTMASKED.fa &&\

bin/kallisto quant -i data/ReferenceGenomes/m_belari_kall_index -o results/kallisto_fragmentlength results/quality_control/trimming/2017_08_08_MRDR5_trim_R1.fastq.gz results/quality_control/trimming/2017_08_08_MRDR5_trim_R2.fastq.gz results/quality_control/trimming/2017_08_08_MRDR6_trim_R1.fastq.gz results/quality_control/trimming/2017_08_08_MRDR6_trim_R2.fastq.gz
