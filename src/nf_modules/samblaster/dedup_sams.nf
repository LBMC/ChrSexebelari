params.sam = "$baseDir/data/sam/*.sam"

log.info "sams files : ${params.sam}"

Channel
  .fromPath( params.sam )
  .ifEmpty { error "Cannot find any sam files matching: ${params.sam}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { sam_files }

process dedup_sam {
  tag "$file_id"
  cpus 4

  input:
    set file_id, file(sam) from sam_files

  output:
    set file_id, "*_dedup.sam*" into dedup_sam_files
  script:
"""
samblaster --addMateTags -i ${sam} -o ${file_id}_dedup.sam
"""
}


