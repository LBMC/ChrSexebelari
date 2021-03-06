./nextflow src/nf_modules/BWA/indexing.nf \
  -c src/nf_modules/BWA/indexing.config \
  -profile docker \
  --fasta "data/tiny_dataset/fasta/tiny_v2.fasta"

# ./nextflow src/nf_modules/BWA/mapping_single.nf \
#   -c src/nf_modules/BWA/mapping_single.config \
#   -profile docker \
#   --index "results/mapping/index/tiny_v2.index" \
#   --fastq "data/tiny_dataset/fastq/tiny*_S.fastq"

./nextflow src/nf_modules/BWA/mapping_paired.nf \
  -c src/nf_modules/BWA/mapping_paired.config \
  -profile docker \
  --index "results/mapping/index/tiny_v2*" \
  --fastq "data/tiny_dataset/fastq/tiny*_R{1,2}.fastq"

