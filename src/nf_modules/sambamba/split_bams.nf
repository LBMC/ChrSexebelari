params.bam = "$baseDir/data/bam/*.bam"

log.info "bams files : ${params.bam}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any bam files matching: ${params.bam}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { bam_files }

process split_bam {
  tag "$file_id"
  cpus 4

  input:
    set file_id, file(bam) from bam_files

  output:
    set file_id, "*_forward.bam*" into forward_bam_files
    set file_id, "*_reverse.bam*" into reverse_bam_files
  script:
"""
sambamba view -t ${task.cpus} -h -F "strand == '+'" ${bam} > ${file_id}_forward.bam
sambamba view -t ${task.cpus} -h -F "strand == '-'" ${bam} > ${file_id}_reverse.bam
"""
}

