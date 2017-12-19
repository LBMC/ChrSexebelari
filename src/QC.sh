# done on local computer 

bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files "data/fastq/*6*" --cpu 4 --fastqc "/usr/local/bin/FastQC/fastqc" --paired true --trimmer "cutadapt" --quality_threshold 20 --pigz_version "2.3" --do_adapter_removal false

bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files "data/fastq/*5*" --cpu 4 --fastqc "/usr/local/bin/FastQC/fastqc" --paired true --trimmer "cutadapt" --urqt "/usr/local/bin/UrQt/UrQt" --quality_threshold 20 --pigz_version "2.3" --do_adapter_removal false
