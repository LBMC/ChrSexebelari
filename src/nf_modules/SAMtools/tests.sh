./nextflow src/nf_modules/SAMtools/sort_bams.nf \
  -c src/nf_modules/SAMtools/sort_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam"

./nextflow src/nf_modules/SAMtools/index_bams.nf \
  -c src/nf_modules/SAMtools/index_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.sort.bam"

./nextflow src/nf_modules/SAMtools/split_bams.nf \
  -c src/nf_modules/SAMtools/split_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam"

./nextflow src/nf_modules/SAMtools/filter_bams.nf \
  -c src/nf_modules/SAMtools/filter_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam" \
  --bed "data/tiny_dataset/OLD/2genes.bed"
