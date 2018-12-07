params.fasta = "$baseDir/data/bam/*.fasta"
params.annotation = "$baseDir/data/bam/*.gff3"

log.info "fasta files : ${params.fasta}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .set { fasta_file }
Channel
  .fromPath( params.annotation )
  .ifEmpty { error "Cannot find any annotation files matching: ${params.annotation}" }
  .set { annotation_file }

process index_fasta {
  tag "$fasta.baseName"
  cpus 4
  publishDir "results/mapping/index/", mode: 'copy'

  input:
    file fasta from fasta_file
    file annotation from annotation_file

  output:
    file "*.index*" into index_files

  script:
  def cmd_annotation = "--gff3 ${annotation}"
  if(annotation ==~ /.*\.gtf$/){
    cmd_annotation = "--gtf ${annotation}"
  }
"""
rsem-prepare-reference -p ${task.cpus} --bowtie2 \
--bowtie2-path \$(which bowtie2 | sed 's/bowtie2\$//g') \
${cmd_annotation} ${fasta} ${fasta.baseName}.index > \
${fasta.baseName}_rsem_bowtie2_report.txt
"""
}


