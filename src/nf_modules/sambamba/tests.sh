./nextflow src/nf_modules/sambamba/sort_bams.nf \
  -c src/nf_modules/sambamba/sort_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam"

./nextflow src/nf_modules/sambamba/index_bams.nf \
  -c src/nf_modules/sambamba/index_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.sort.bam"

./nextflow src/nf_modules/sambamba/split_bams.nf \
  -c src/nf_modules/sambamba/split_bams.config \
  -profile docker \
  --bam "data/tiny_dataset/map/tiny_v2.bam"
