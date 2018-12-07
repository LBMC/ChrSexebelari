./nextflow src/nf_modules/HTSeq/htseq.nf \
  -c src/nf_modules/HTSeq/htseq.config \
  -profile docker \
  --gtf "data/tiny_dataset/annot/tiny.gff" \
  --bam "data/tiny_dataset/map/tiny_v2.bam"

