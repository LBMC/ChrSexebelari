params.fasta = "$baseDir/data/bam/*.fasta"
params.fastq = "$baseDir/data/fastq/*_{1,2}.fastq"
params.contig = ""

log.info "fasta files : ${params.fasta}"
log.info "fastq files : ${params.fastq}"
log.info "contig files : ${params.contig}"

def contig_list = Eval.me(params.contig)

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any bam files matching: ${params.fasta}" }
  .splitFasta( record: [id: true, text: true])
  .filter { record -> record.id in contig_list }
  .map{record -> record.text}
  .collectFile( name:'sequences.fasta')
  .set { contig_file }

Channel
  .fromFilePairs( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process contigs_fasta {
  publishDir "results/mapping/fasta/", mode: 'copy'
  input:
    file fasta from contig_file
  output:
    file "*.fasta" into fasta_file

  script:
"""
cp ${fasta} contigs.fasta
"""
}

process index_fasta {
  tag "$fasta.baseName"
  cpus 4
  publishDir "results/mapping/index/", mode: 'copy'

  input:
    file fasta from fasta_file

  output:
    file "*.index*" into index_files
    file "*_report.txt" into indexing_report

  script:
"""
bowtie2-build --threads ${task.cpus} ${fasta} ${fasta.baseName}.index &> ${fasta.baseName}_bowtie2_report.txt

if grep -q "Error" ${fasta.baseName}_bowtie2_report.txt; then
  exit 1
fi
"""
}

process mapping_fastq {
  tag "$pair_id"
  cpus 4

  input:
  set pair_id, file(reads) from fastq_files
  file index from index_files.collect()

  output:
  set pair_id, "*.bam" into bam_files
  file "*_report.txt" into mapping_report

  script:
  index_id = index[0]
  for (index_file in index) {
    if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
        index_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
    }
  }
"""
bowtie2 --very-sensitive -p ${task.cpus} -x ${index_id} \
-1 ${reads[0]} -2 ${reads[1]} 2> \
${pair_id}_bowtie2_report.txt | \
samtools view -Sb - > ${pair_id}.bam

if grep -q "Error" ${pair_id}_bowtie2_report.txt; then
  exit 1
fi
"""
}

process sort_bam {
  tag "$file_id"
  cpus 4
  publishDir "results/mapping/bams/", mode: 'copy'

  input:
    set file_id, file(bam) from bam_files

  output:
    set file_id, "*_sorted.bam" into sorted_bam_files

  script:
"""
sambamba sort -t ${task.cpus} -o ${file_id}_sorted.bam ${bam}
"""
}

process index_bam {
  tag "$file_id"
  cpus 4
  publishDir "results/mapping/bams/", mode: 'copy'

  input:
    set file_id, file(bam) from sorted_bam_files

  output:
    set file_id, "*.bam*" into indexed_bam_file

  script:
"""
sambamba index -t ${task.cpus} ${bam}
"""
}

