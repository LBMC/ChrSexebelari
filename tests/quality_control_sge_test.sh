#!/bin/sh

module load nextflow/0.25.1

nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'cutadapt' --cpu 4
nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --trimmer 'cutadapt' --cpu 4  && \
nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'UrQt' --cpu 4  && \
nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --trimmer 'UrQt' --cpu 4

ls -l data/examples/tiny_dataset/fastq/*.fastq | awk '{system("gzip "$9" -c > "$9".gz")}'

nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq.gz' --trimmer 'cutadapt' --cpu 4
nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq.gz' --trimmer 'cutadapt' --cpu 4  && \
nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq.gz' --trimmer 'UrQt' --cpu 4  && \
nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq.gz' --trimmer 'UrQt' --cpu 4
