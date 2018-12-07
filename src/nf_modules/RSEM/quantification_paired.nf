params.fastq = "$baseDir/data/fastq/*_{1,2}.fastq"
params.index = "$baseDir/data/index/*.index.*"

log.info "fastq files : ${params.fastq}"
log.info "index files : ${params.index}"

Channel
  .fromFilePairs( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }
Channel
  .fromPath( params.index )
  .ifEmpty { error "Cannot find any index files matching: ${params.index}" }
  .set { index_files }

process mapping_fastq {
  tag "$pair_id"
  cpus 4
  publishDir "results/mapping/quantification/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files
  file index from index_files.toList()

  output:
  file "*" into counts_files

  script:
  index_id = index[0]
  for (index_file in index) {
    if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
        index_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
    }
  }
"""
rsem-calculate-expression --bowtie2 \
--bowtie2-path \$(which bowtie2 | sed 's/bowtie2\$//g') \
--bowtie2-sensitivity-level "very_sensitive" \
-output-genome-bam -p ${task.cpus} \
--paired-end ${reads[0]} ${reads[1]} ${index_id} ${pair_id} \
2> ${pair_id}_rsem_bowtie2_report.txt

if grep -q "Error" ${pair_id}_rsem_bowtie2_report.txt; then
  exit 1
fi
"""
}


