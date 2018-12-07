./nextflow src/nf_modules/cutadapt/adaptor_removal_paired.nf \
  -c src/nf_modules/cutadapt/adaptor_removal_paired.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq"

./nextflow src/nf_modules/cutadapt/adaptor_removal_single.nf \
  -c src/nf_modules/cutadapt/adaptor_removal_single.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq"

./nextflow src/nf_modules/cutadapt/trimming_paired.nf \
  -c src/nf_modules/cutadapt/trimming_paired.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq"

./nextflow src/nf_modules/cutadapt/trimming_single.nf \
  -c src/nf_modules/cutadapt/trimming_single.config \
  -profile docker \
  --fastq "data/tiny_dataset/fastq/tiny*_S.fastq"
