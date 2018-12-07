./nextflow src/nf_modules/MultiQC/multiqc_paired.nf \
  -c src/nf_modules/MultiQC/multiqc_paired.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq"

./nextflow src/nf_modules/MultiQC/multiqc_single.nf \
  -c src/nf_modules/MultiQC/multiqc_single.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_S.fastq"
