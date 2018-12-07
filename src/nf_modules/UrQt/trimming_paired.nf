log.info "fastq files : ${params.fastq}"

Channel
  .fromFilePairs( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process trimming {
  tag "${reads}"
  cpus 4
  publishDir "results/fastq/trimming/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files

  output:
  set pair_id, "*_trim_R{1,2}.fastq.gz" into fastq_files_trim

  script:
"""
UrQt --t 20 --m ${task.cpus} --gz \
--in ${reads[0]} --inpair ${reads[1]} \
--out ${pair_id}_trim_R1.fastq.gz --outpair ${pair_id}_trim_R2.fastq.gz \
> ${pair_id}_trimming_report.txt
"""
}

