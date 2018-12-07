./nextflow src/nf_modules/samblaster/dedup_sams.nf \
  -c src/nf_modules/samblaster/dedup_sams.config \
  -profile docker \
  --sam "data/tiny_dataset/map/tiny_v2.sam"
