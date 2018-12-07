params.bam = "$baseDir/data/bam/*.bam"

log.info "bams files : ${params.bam}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any bam files matching: ${params.bam}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { bam_files }

process index_bam {
  tag "$file_id"
  cpus 4

  input:
    set file_id, file(bam) from bam_files

  output:
    set file_id, "*.bam*" into indexed_bam_file

  script:
"""
sambamba index -t ${task.cpus} ${bam}
"""
}

