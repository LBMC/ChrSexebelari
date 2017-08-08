# This will call quality control profile in nextflow.config in local mode.

#RootDir is: 
#{/home/cburny/Bureau/Projects/ChrSexe/src/RNASeq/}


sudo bin/nextflow src/RNASeq/src/pipe/quality_control.nf -c src/RNASeq/src/nextflow.config -profile quality_control --fastq_files "/home/cburny/Bureau/Projects/ChrSexe/data/PoolSeq/fastq/*MRDR5_R{1,2}.fastq.gz" --cpu 2 --fastqc "/usr/local/bin/FastQC/fastqc" --trimmer "cutadapt" --urqt "/usr/local/bin/UrQt/UrQt" --quality_threshold 20 --pigz_version "2.3" --do_adapter_removal false --pref_multiqc "MRDR5"

sudo bin/nextflow src/RNASeq/src/pipe/quality_control.nf -c src/RNASeq/src/nextflow.config -profile quality_control --fastq_files "/home/cburny/Bureau/Projects/ChrSexe/data/PoolSeq/fastq/*MRDR6_R{1,2}.fastq.gz" --cpu 2 --fastqc "/usr/local/bin/FastQC/fastqc" --trimmer "cutadapt" --urqt "/usr/local/bin/UrQt/UrQt" --quality_threshold 20 --pigz_version "2.3" --do_adapter_removal false --pref_multiqc "MRDR6"

