params.bam = "$baseDir/data/bam/*.bam"

log.info "bams files : ${params.bam}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any bam files matching: ${params.bam}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { bam_files }

process sort_bam {
  tag "$file_id"
  cpus 4

  input:
    set file_id, file(bam) from bam_files

  output:
    set file_id, "*_sorted.bam" into sorted_bam_files

  script:
"""
samtools sort -@ ${task.cpus} -O BAM -o ${file_id}_sorted.bam ${bam}
"""
}

