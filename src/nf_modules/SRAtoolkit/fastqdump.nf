/*
* sra-tools :

*/

/*                      fastq-dump
* Imputs : srr list
* Outputs : fastq files
*/

params.list_srr = "$baseDir/data/SRR/*.txt"

log.info "downloading list srr : ${params.list_srr}"

Channel
  .fromPath( params.list_srr )
  .ifEmpty { error "Cannot find any bam files matching: ${params.list_srr}" }
  .splitCsv()
  .map { it -> it[0]}
  .set { SRR }

//run is the column name containing SRR ids

process fastq_dump {
  tag "$file_id"
  publishDir "results/download/fastq/${file_id}/", mode: 'copy'

  input:
    val file_id from SRR

  output:
    set file_id, "*.fastq" into fastq

  script:
"""
#for test only 10000  reads are downloading with the option -N 10000 -X 20000
fastq-dump --split-files --defline-seq '@\$ac_\$si/\$ri' --defline-qual "+" -N 10000 -X 20000 ${file_id}
if [ -f ${file_id}_1.fastq ]
then
  mv ${file_id}_1.fastq ${file_id}_R1.fastq
fi
if [ -f ${file_id}_2.fastq ]
then
  mv ${file_id}_2.fastq ${file_id}_R2.fastq
fi
"""
}
