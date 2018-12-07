params.bam = "$baseDir/data/bam/*.bam"
params.gtf = "$baseDir/data/annotation/*.gtf"

log.info "bam files : ${params.bam}"
log.info "gtf files : ${params.gtf}"

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.bam}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { bam_files }
Channel
  .fromPath( params.gtf )
  .ifEmpty { error "Cannot find any gtf file matching: ${params.gtf}" }
  .set { gtf_file }

process sort_bam {
  tag "$file_id"
  cpus 4

  input:
    set file_id, file(bam) from bam_files

  output:
    set file_id, "*_sorted.sam" into sorted_bam_files

  script:
"""
# sort bam by name
samtools sort -@ ${task.cpus} -n -O SAM -o ${file_id}_sorted.sam ${bam}
"""
}

process counting {
  tag "$file_id"
  publishDir "results/quantification/", mode: 'copy'

  input:
  set file_id, file(bam) from sorted_bam_files
  file gtf from gtf_file

  output:
  file "*.count" into count_files

  script:
"""
htseq-count ${bam} ${gtf} \
-r pos --mode=intersection-nonempty -a 10 -s no -t exon -i gene_id \
> ${file_id}.count
"""
}

