params.bam = "$baseDir/data/bam/*.bam"
params.bed = "$baseDir/data/bam/*.bed"

log.info "bams files : ${params.bam}"
log.info "bed file : ${params.bed}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any bam files matching: ${params.bam}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { bam_files }
Channel
  .fromPath( params.bed )
  .ifEmpty { error "Cannot find any bed file matching: ${params.bed}" }
  .set { bed_files }

process filter_bam {
  tag "$file_id"
  cpus 4

  input:
    set file_id, file(bam) from bam_files
    file bed from bed_files

  output:
    set file_id, "*_filtered.bam*" into filtered_bam_files
  script:
"""
samtools view -@ ${task.cpus} -hb ${bam} -L ${bed} > ${file_id}_filtered.bam
"""
}


