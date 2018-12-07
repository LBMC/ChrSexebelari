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
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .groupTuple()
  .set { index_files }

process mapping_fastq {
  tag "$reads"
  cpus 4
  publishDir "results/mapping/sam/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files
  set index_id, file(index) from index_files.collect()

  output:
  file "${pair_id}.sam" into sam_files
  file "${pair_id}_bwa_report.txt" into mapping_repport_files

  script:
"""
bwa mem -t ${task.cpus} \
${index_id} ${reads[0]} ${reads[1]} \
-o ${pair_id}.sam &> ${pair_id}_bwa_report.txt
"""
}

