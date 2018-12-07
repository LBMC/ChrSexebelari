cp results/training/bams/sBNLN18.bam results/training/bams/sBNLN18_control.bam
./nextflow src/nf_modules/MUSIC/peak_calling_single.nf \
  -c src/nf_modules/MUSIC/peak_calling_single.config \
  -profile docker \
  --fasta "results/training/fasta/*.fasta" \
  --bam "results/training/bams/s*.bam" \
  --index "results/training/mapping/index/*" \
  --read_size 50 --frag_size 300
