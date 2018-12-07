params.fasta = "$baseDir/data/bam/*.fasta"

log.info "fasta files : ${params.fasta}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any bam files matching: ${params.fasta}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { fasta_file }

process index_fasta {
  tag "$fasta_id"
  cpus 4
  publishDir "results/mapping/index/", mode: 'copy'

  input:
    set fasta_id, file(fasta) from fasta_file

  output:
    set fasta_id, "${fasta.baseName}.*" into index_files
    file "*_bwa_report.txt" into index_files_report

  script:
"""
bwa index -p ${fasta.baseName} ${fasta} \
&> ${fasta.baseName}_bwa_report.txt
"""
}

