params.fastq = "$baseDir/data/fastq/*.fastq"
params.index = "$baseDir/data/index/*.index*"
params.mean = 200
params.sd = 100

log.info "fastq files : ${params.fastq}"
log.info "index files : ${params.index}"
log.info "mean read size: ${params.mean}"
log.info "sd read size: ${params.sd}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { fastq_files }
Channel
  .fromPath( params.index )
  .ifEmpty { error "Cannot find any index files matching: ${params.index}" }
  .set { index_files }

process mapping_fastq {
  tag "$file_id"
  cpus 4
  publishDir "results/mapping/quantification/", mode: 'copy'

  input:
  set file_id, file(reads) from fastq_files
  file index from index_files.toList()

  output:
  file "*" into count_files

  script:
  index_id = index[0]
  for (index_file in index) {
    if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
        index_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
    }
  }
"""
ls -l
rsem-calculate-expression --bowtie2 \
--bowtie2-path \$(which bowtie2 | sed 's/bowtie2\$//g') \
--bowtie2-sensitivity-level "very_sensitive" \
--fragment-length-mean ${params.mean} --fragment-length-sd ${params.sd} \
--output-genome-bam -p ${task.cpus} \
${reads} ${index_id} ${file_id} \
2> ${file_id}_rsem_bowtie2_report.txt

if grep -q "Error" ${file_id}_rsem_bowtie2_report.txt; then
  exit 1
fi
"""
}

