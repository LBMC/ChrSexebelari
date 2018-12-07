params.fastq = "$baseDir/data/*.fastq"
params.fasta = "$baseDir/data/*.fasta"
params.seq_length = 800000
log.info "fastq files : ${params.fastq}"
log.info "fasta files : ${params.fasta}"
log.info "fasta length to retain : ${params.seq_length}"
def normal_sample = Eval.me(params.normal)
def tumor_sample = Eval.me(params.tumor)
log.info "normal : ${normal_sample}"
log.info "tumor : ${tumor_sample}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { fasta_file  }
Channel
  .fromFilePairs( params.fastq )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process adaptor_removal {
  tag "$pair_id"
  publishDir "results/fastq/adaptor_removal/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files

  output:
  set pair_id, file("*.fastq.gz") into fastq_files_cut
  file "*_cutadapt_report.txt" into cut_files_report

  script:
"""
cutadapt -a AGATCGGAAGAG -g CTCTTCCGATCT -A AGATCGGAAGAG -G CTCTTCCGATCT \
-o ${pair_id}_cut_R1.fastq.gz -p ${pair_id}_cut_R2.fastq.gz \
${reads[0]} ${reads[1]} > ${pair_id}_cutadapt_report.txt
"""
}

process trimming {
  tag "${pair_id}"
  cpus 4
  publishDir "results/fastq/trimming/", mode: 'copy'

  input:
  set pair_id, file(reads) from fastq_files_cut

  output:
  set pair_id, file("*.fastq.gz") into fastq_files_trim
  file "*_trimming_report.txt" into trimming_files_report

  script:
"""
UrQt --t 20 --m ${task.cpus} --gz \
--in ${reads[0]} --inpair ${reads[1]} \
--out ${pair_id}_trim_R1.fastq.gz --outpair ${pair_id}_trim_R2.fastq.gz \
> ${pair_id}_trimming_report.txt
"""
}

process filter_fasta {
  tag "$fasta_id"
  cpus 4
  publishDir "results/fasta/", mode: 'copy'

  input:
    set fasta_id, file(fasta) from fasta_file

  output:
    set fasta_idf, "*_filtered.fasta" into filter_fasta_files

  script:
    fasta_idf = "${fasta_id}_filtered"
"""
bioawk -c fastx '{ if(length(\$seq) > $params.seq_length) { print ">"\$name; print \$seq }}' ${fasta} > \
${fasta_id}_filtered.fasta
"""
}

filter_fasta_files.into{
  filtered_fasta_files;
  indel_fasta_file;
  recalibration_fasta_file;
  haplotypecaller_fasta_file
}

process index_fasta {
  tag "$file_id"
  cpus 12
  publishDir "results/mapping/index/", mode: 'copy'

  input:
    set file_id, file(fasta) from filtered_fasta_files

  output:
    file "*.index*" into index_files
    file "*_report.txt" into indexing_report

  script:
"""
bowtie2-build --threads ${task.cpus} ${fasta} ${file_id}.index &> ${file_id}_bowtie2_report.txt

if grep -q "Error" ${file_id}_bowtie2_report.txt; then
  exit 1
fi
"""
}

fastq_files_trim.into{
  fastq_files_trim_norm;
  fastq_files_trim_tumor
}

collect_fastq_files_trim_norm = fastq_files_trim_norm
  .filter{ normal_sample.contains(it[0]) }
  .map { it -> ["normal_sample", it[0], it[1]]}

collect_fastq_files_trim_tumor = fastq_files_trim_tumor
  .filter{ tumor_sample.contains(it[0]) }
  .map { it -> ["tumor_sample", it[0], it[1]]}

collect_fastq_files_trim = Channel.create()
  .mix(collect_fastq_files_trim_norm, collect_fastq_files_trim_tumor)

process mapping_fastq {
  tag "$pair_id"
  cpus 12
  publishDir "results/mapping/bam/", mode: 'copy'

  input:
  set sample_name, pair_id, file(reads) from collect_fastq_files_trim
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
--rg-id ${sample_name} \
--rg PL:Illumina \
--rg SM:${sample_name} \
-1 ${reads[0]} -2 ${reads[1]} 2> \
${pair_id}_bowtie2_report.txt | \
samblaster --addMateTags -M -i /dev/stdin | \
sambamba view -t ${task.cpus} --valid -S -f bam -l 0 /dev/stdin \
-o ${pair_id}.bam

if grep -q "Error" ${pair_id}_bowtie2_report.txt; then
  exit 1
fi
"""
}

process sort_bam {
  tag "$file_id"
  cpus 4

  input:
    set file_id, file(bam) from bam_files

  output:
    set file_id, "*_sorted.bam" into sorted_bam_files

  script:
"""
sambamba sort -t ${task.cpus} --tmpdir=./tmp -o ${file_id}_sorted.bam ${bam}
"""
}

sorted_bam_files.into {
  sorted_bam_file_norm;
  sorted_bam_file_tumor;
}

collect_sorted_bam_file_norm = sorted_bam_file_norm
  .filter{ normal_sample.contains(it[0]) }
  .map { it -> it[1]}
  .buffer( size: normal_sample.size())
  .map { it -> ["normal_sample", it]}
collect_sorted_bam_file_tumor = sorted_bam_file_tumor
  .filter{ tumor_sample.contains(it[0]) }
  .map { it -> it[1]}
  .buffer( size: tumor_sample.size())
  .map { it -> ["tumor_sample", it]}

collect_sorted_bam_file = Channel.create()
  .mix(collect_sorted_bam_file_norm, collect_sorted_bam_file_tumor)

process merge_bam {
  tag "$file_id"
  cpus 4
  publishDir "results/mapping/bam/", mode: 'copy'

  input:
    set file_id, file(bam) from collect_sorted_bam_file

  output:
    set file_id, "*.bam" into merged_bam_files

  script:
"""
if ((\$(ls -l *.bam | wc -l) > 1)); then
sambamba merge -t ${task.cpus} ${file_id}.bam ${bam}
else
cp ${bam} ${file_id}.bam
fi
"""
}

merged_bam_files.into{
  index_merged_bam_files;
  haplo_bam_files_norm;
  haplo_bam_files_tumor
}

process index_bam {
  tag "$file_id"
  cpus 4
  publishDir "results/mapping/bam/", mode: 'copy'

  input:
    set file_id, file(bam) from index_merged_bam_files

  output:
    set file_id, "*.bam.bai" into index_bam_files

  script:
"""
sambamba index -t ${task.cpus} ${bam}
"""
}

index_bam_files.into{
  named_index_bam_files;
  indexed_bam_files
}

haplotypecaller_fasta_file.into{
    final_fasta_file;
    index2_fasta_file
    index3_fasta_file
  }

process index2_fasta {
  tag "$genome_id"
  publishDir "results/fasta/", mode: 'copy'

  input:
    set genome_id, file(fasta) from index2_fasta_file

  output:
    set genome_id, "*.dict" into indexed2_fasta_file

  script:
"""
gatk CreateSequenceDictionary -R ${fasta} &> gatk_output.txt
"""
}

process index3_fasta {
  tag "$genome_id"
  publishDir "results/fasta/", mode: 'copy'

  input:
    set genome_id, file(fasta) from index3_fasta_file

  output:
    set genome_id, "*.fai" into indexed3_fasta_file

  script:
"""
samtools faidx ${fasta}
"""
}

final_bam_files_norm = haplo_bam_files_norm
  .filter{ "normal_sample" == it[0] }
final_bam_files_tumor = haplo_bam_files_tumor
  .filter{ "tumor_sample" == it[0] }

indexed_bam_files.into {
  index_bam_files_norm;
  index_bam_files_tumor
}
final_indexed_bam_files_norm = index_bam_files_norm
  .filter{ "normal_sample" == it[0] }
final_indexed_bam_files_tumor = index_bam_files_tumor
   .filter{ "tumor_sample" == it[0] }

final_bam_files_norm.set{
  samtools_SNP_bam_files_norm
}
final_bam_files_tumor.set{
  samtools_SNP_bam_files_tumor;
}
final_indexed_bam_files_norm.set{
  samtools_SNP_index_bam_files_norm
}
final_indexed_bam_files_tumor.set{
  samtools_SNP_index_bam_files_tumor;
}
final_fasta_file.into{
  samtools_SNP_fasta_file_tumor;
  samtools_SNP_fasta_file_norm;
}
indexed2_fasta_file.into{
  samtools_SNP_indexed2_fasta_file_tumor;
  samtools_SNP_indexed2_fasta_file_norm;
}
indexed3_fasta_file.into{
  samtools_SNP_indexed3_fasta_file_tumor;
  samtools_SNP_indexed3_fasta_file_norm;
}

process samtools_SNP_tumor {
  tag "$file_id_tumor"
  cpus 1
  publishDir "results/SNP/vcf_samtools/", mode: 'copy'

  input:
    set file_id_tumor, file(bam_tumor) from samtools_SNP_bam_files_tumor
    set file_ididx_tumor, file(bamidx_tumor) from samtools_SNP_index_bam_files_tumor
    set genome_id, file(fasta) from samtools_SNP_fasta_file_tumor
    set genome2_idx, file(fasta2idx) from samtools_SNP_indexed2_fasta_file_tumor
    set genome3_idx, file(fasta3idx) from samtools_SNP_indexed3_fasta_file_tumor

  output:
    set file_id_tumor, "*.vcf" into vcf_files_tumor

  script:
"""
bcftools mpileup -AE -f ${fasta} ${bam_tumor} --output-type v \
-a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR | \
bcftools call -mv --output-type v > ${file_id_tumor}_raw.vcf
bcftools filter -s LowQual -e '%QUAL<20 || DP>100' ${file_id_tumor}_raw.vcf \
> ${file_id_tumor}_filtered.vcf
"""
}

process samtools_SNP_norm {
  tag "$file_id_norm"
  cpus 1
  publishDir "results/SNP/vcf_samtools/", mode: 'copy'

  input:
    set file_id_norm, file(bam_norm) from samtools_SNP_bam_files_norm
    set file_ididx_norm, file(bamidx_norm) from samtools_SNP_index_bam_files_norm
    set genome_id, file(fasta) from samtools_SNP_fasta_file_norm
    set genome2_idx, file(fasta2idx) from samtools_SNP_indexed2_fasta_file_norm
    set genome3_idx, file(fasta3idx) from samtools_SNP_indexed3_fasta_file_norm

  output:
    set file_id_norm, "*.vcf" into vcf_files_norm

  script:
"""
bcftools mpileup -AE -f ${fasta} ${bam_norm} --output-type v \
-a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR | \
bcftools call -mv --output-type v  > ${file_id_norm}_raw.vcf
bcftools filter -s LowQual -e '%QUAL<20 || DP>100' ${file_id_norm}_raw.vcf \
> ${file_id_norm}_filtered.vcf
"""
}

process vcf_to_csv_tumor {
  tag "$file_id_tumor"
  cpus 1
  publishDir "results/SNP/vcf_samtools/", mode: 'copy'

  input:
    set file_id_tumor, file(vcf) from vcf_files_tumor

  output:
    set file_id_tumor, "*.csv" into csv_files_tumor

  script:
"""
gatk VariantsToTable -V ${file_id_tumor}_raw.vcf \
-F CHROM -F POS -F TYPE -GF GT -GF AD -F AD -F DP \
-O ${file_id_tumor}_raw.csv
gatk VariantsToTable -V ${file_id_tumor}_filtered.vcf \
-F CHROM -F POS -F TYPE -GF GT -GF AD -F AD -F DP \
-O ${file_id_tumor}_filtered.csv
"""
}

process vcf_to_csv_norm {
  tag "$file_id_norm"
  cpus 1
  publishDir "results/SNP/vcf_samtools/", mode: 'copy'

  input:
    set file_id_norm, file(vcf) from vcf_files_norm

  output:
    set file_id_norm, "*.csv" into csv_files_norm

  script:
"""
gatk VariantsToTable -V ${file_id_norm}_raw.vcf \
-F CHROM -F POS -F TYPE -GF GT -GF AD -F AD -F DP \
-O ${file_id_norm}_raw.csv
gatk VariantsToTable -V ${file_id_norm}_filtered.vcf \
-F CHROM -F POS -F TYPE -GF GT -GF AD -F AD -F DP \
-O ${file_id_norm}_filtered.csv
"""
}

