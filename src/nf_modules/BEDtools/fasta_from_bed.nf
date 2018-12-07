/*
* bedtools :
* Imputs : fasta files
* Imputs : bed files
* Output : fasta files
*/
/*                      fasta extraction                                     */

params.fasta = "$baseDir/data/fasta/*.fasta"
params.bed = "$baseDir/data/annot/*.bed"

log.info "fasta file : ${params.fasta}"
log.info "bed file : ${params.bed}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .set { fasta_files }
Channel
  .fromPath( params.bed )
  .ifEmpty { error "Cannot find any bed files matching: ${params.bed}" }
  .set { bed_files }

process fasta_from_bed {
  tag "${bed.baseName}"
  cpus 4
  publishDir "results/fasta/", mode: 'copy'

  input:
  file fasta from fasta_files
  file bed from bed_files

  output:
  file "*_extracted.fasta" into fasta_files_extracted

  script:
"""
bedtools getfasta -name \
-fi ${fasta} -bed ${bed} -fo ${bed.baseName}_extracted.fasta
"""
}
