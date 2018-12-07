params.fastq = "$baseDir/data/fastq/*.fastq"

log.info "fastq files : ${params.fastq}"

Channel
  .fromPath( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { fastq_files }

process fastqc_fastq {
  tag "$file_id"
  publishDir "results/fastq/fastqc/", mode: 'copy'
  cpus = 1

  input:
    set file_id, file(reads) from fastq_files

  output:
    file "*.{zip,html}" into fastqc_report

  script:
"""
fastqc --quiet --threads ${task.cpus} --format fastq --outdir ./ ${reads}
"""
}

process multiqc {
  tag "$report[0].baseName"
  publishDir "results/fastq/multiqc/", mode: 'copy'
  cpus = 1

  input:
    file report from fastqc_report.collect()

  output:
    file "*multiqc_*" into multiqc_report

  script:
"""
multiqc -f .
"""
}

