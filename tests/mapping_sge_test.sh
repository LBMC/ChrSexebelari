#!/bin/sh
module load nextflow/0.25.1

cp data/examples/tiny_dataset/fastq/tiny_R1.fastq data/examples/tiny_dataset/fastq/c1tiny_R1.fastq
cp data/examples/tiny_dataset/fastq/tiny_R2.fastq data/examples/tiny_dataset/fastq/c1tiny_R2.fastq
cp data/examples/tiny_dataset/fastq/tiny_R1.fastq data/examples/tiny_dataset/fastq/c2tiny_R1.fastq
cp data/examples/tiny_dataset/fastq/tiny_R2.fastq data/examples/tiny_dataset/fastq/c2tiny_R2.fastq

nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --mapper 'salmon' --paired true --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff' --cpu 16 && \
nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --mapper 'salmon' --paired false --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 16 && \
nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --mapper 'kallisto' --paired true --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 16 && \
nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --mapper 'kallisto' --paired false --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 16 && \
nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --mapper 'bowtie2' --paired true --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 16 && \
nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --mapper 'bowtie2' --paired false --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 16 && \
nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --mapper 'bowtie2' --paired true --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gtf' --quantifier 'rsem' --cpu 16


ls -l data/examples/tiny_dataset/fastq/*.fastq | awk '{system("gzip "$9" -c > "$9".gz")}'

nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq.gz' --mapper 'salmon' --paired true --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff' --cpu 16 && \
nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq.gz' --mapper 'salmon' --paired false --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 16 && \
nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq.gz' --mapper 'kallisto' --paired true --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 16 && \
nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq.gz' --mapper 'kallisto' --paired false --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 16 && \
nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq.gz' --mapper 'bowtie2' --paired true --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 16 && \
nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq.gz' --mapper 'bowtie2' --paired false --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 16 && \
nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_sge --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq.gz' --mapper 'bowtie2' --paired true --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gtf' --quantifier 'rsem' --cpu 16
