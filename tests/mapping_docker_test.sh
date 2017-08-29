#!/bin/sh
cp data/examples/tiny_dataset/fastq/tiny_R1.fastq data/examples/tiny_dataset/fastq/c1tiny_R1.fastq
cp data/examples/tiny_dataset/fastq/tiny_R2.fastq data/examples/tiny_dataset/fastq/c1tiny_R2.fastq
cp data/examples/tiny_dataset/fastq/tiny_R1.fastq data/examples/tiny_dataset/fastq/c2tiny_R1.fastq
cp data/examples/tiny_dataset/fastq/tiny_R2.fastq data/examples/tiny_dataset/fastq/c2tiny_R2.fastq

cp src/file_handle/src/file_handle.py src/pipe/mapping/ && \
docker build src/pipe/mapping -t 'mapping:0.0.1' && \

bin/nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --mapper 'salmon' --paired true --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff' --cpu 2 && \
bin/nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --mapper 'salmon' --paired false --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 2 && \
bin/nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --mapper 'kallisto' --paired true --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 2 && \
bin/nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --mapper 'kallisto' --paired false --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 2 && \
bin/nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --mapper 'bowtie2' --paired true --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 2 && \
bin/nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --mapper 'bowtie2' --paired false --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 2 && \
bin/nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --mapper 'bowtie2' --paired true --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gtf' --quantifier 'rsem' --cpu 2

ls -l data/examples/tiny_dataset/fastq/*.fastq | awk '{system("gzip "$9" -c > "$9".gz")}'

bin/nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq.gz' --mapper 'salmon' --paired true --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff' --cpu 2 && \
bin/nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq.gz' --mapper 'salmon' --paired false --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 2 && \
bin/nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq.gz' --mapper 'kallisto' --paired true --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 2 && \
bin/nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq.gz' --mapper 'kallisto' --paired false --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 2 && \
bin/nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq.gz' --mapper 'bowtie2' --paired true --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 2 && \
bin/nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq.gz' --mapper 'bowtie2' --paired false --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gff'  --cpu 2 && \
bin/nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq.gz' --mapper 'bowtie2' --paired true --reference 'data/examples/tiny_dataset/fasta/*tiny_v2.fasta' --annotation 'data/examples/tiny_dataset/annot/*tiny.gtf' --quantifier 'rsem' --cpu 2
