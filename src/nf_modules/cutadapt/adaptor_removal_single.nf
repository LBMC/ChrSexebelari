log.info "fastq files : ${params.fastq}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { fastq_files }

process adaptor_removal {
  tag "$file_id"

  input:
  set file_id, file(reads) from fastq_files

  output:
  set file_id, "*_cut.fastq.gz" into fastq_files_cut

  script:
  """
  cutadapt -a AGATCGGAAGAG -g CTCTTCCGATCT\
  -o ${file_id}_cut.fastq.gz \
  ${reads} > ${file_id}_report.txt
  """
}

