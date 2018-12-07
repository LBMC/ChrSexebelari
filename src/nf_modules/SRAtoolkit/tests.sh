./nextflow src/nf_modules/SRAtoolkit/fastqdump.nf \
  -c src/nf_modules/SRAtoolkit/fastqdump.config \
  -profile docker \
  --list_srr "src/nf_modules/SRAtoolkit/list-srr.txt"
