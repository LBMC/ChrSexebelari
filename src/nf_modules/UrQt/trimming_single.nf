params.fastq = "$baseDir/data/fastq/*.fastq"

log.info "fastq files : ${params.fastq}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { fastq_files }

process trimming {
  tag "$file_id"
  cpus 4

  input:
  set file_id, file(reads) from fastq_files

  output:
  set file_id, "*_trim.fastq.gz" into fastq_files_trim

  script:
  """
  UrQt --t 20 --m ${task.cpus} --gz \
  --in ${reads} \
  --out ${file_id}_trim.fastq.gz \
  > ${file_id}_trimming_report.txt
  """
}

