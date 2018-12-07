params.fasta = "$baseDir/data/bam/*.fasta"

log.info "fasta files : ${params.fasta}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any bam files matching: ${params.fasta}" }
  .set { fasta_file }

process index_fasta {
  tag "$fasta.baseName"
  cpus 4
  publishDir "results/mapping/index/", mode: 'copy'

  input:
    file fasta from fasta_file

  output:
    file "*.index*" into index_files
    file "*_kallisto_report.txt" into index_files_report

  script:
"""
kallisto index -k 31 --make-unique -i ${fasta.baseName}.index ${fasta} \
2> ${fasta.baseName}_kallisto_report.txt
"""
}

