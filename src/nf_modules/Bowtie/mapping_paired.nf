/*
* mapping paired fastq
*/

params.fastq = "$baseDir/data/fastq/*_{1,2}.fastq"
params.index = "$baseDir/data/index/*.index.*"

log.info "fastq files : ${params.fastq}"
log.info "index files : ${params.index}"

Channel
  .fromFilePairs( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }
Channel
  .fromPath( params.index )
  .ifEmpty { error "Cannot find any index files matching: ${params.index}" }
  .set { index_files }

process mapping_fastq {
  tag "$pair_id"
  cpus 4
  publishDir "results/mapping/bams/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files
  file index from index_files.collect()

  output:
  file "*.bam" into bam_files
  file "*_report.txt" into mapping_report

  script:
  index_id = index[0]
  for (index_file in index) {
  if (index_file =~ /.*\.1\.ebwt/ && !(index_file =~ /.*\.rev\.1\.ebwt/)) {
        index_id = ( index_file =~ /(.*)\.1\.ebwt/)[0][1]
    }
  }
"""
# -v specify the max number of missmatch, -k the number of match reported per
# reads
bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_id} \
-1 ${reads[0]} -2 ${reads[1]} 2> \
${pair_id}_bowtie_report.txt | \
samtools view -Sb - > ${pair_id}.bam

if grep -q "Error" ${pair_id}_bowtie_report.txt; then
  exit 1
fi
"""
}


