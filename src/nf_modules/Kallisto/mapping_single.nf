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
  file index from index_files.collect()

  output:
  file "*" into count_files

  script:
"""
mkdir ${file_id}
kallisto quant -i ${index} -t ${task.cpus} --single \
--bias --bootstrap-samples 100 -o ${file_id} \
-l ${params.mean} -s ${params.sd} \
${reads} &> ${file_id}/kallisto_report.txt
"""
}

