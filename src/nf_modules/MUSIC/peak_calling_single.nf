params.read_size = 100
params.frag_size = 200
params.step_l = 50
params.min_l = 200
params.max_l = 5000
log.info "bam files : ${params.bam}"
log.info "index files : ${params.index}"
log.info "fasta files : ${params.fasta}"

Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any bam files matching: ${params.fasta}" }
  .set { fasta_files }
Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any bam files matching: ${params.bam}" }
  .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
  .set { bam_files }
Channel
  .fromPath( params.index )
  .ifEmpty { error "Cannot find any index files matching: ${params.index}" }
  .set { index_files }

process compute_mappability {
  tag "${fasta.baseName}"

  input:
    file index from index_files.collect()
    file fasta from fasta_files

  output:
    file "*.bin" into mappability
    file "temp/chr_ids.txt" into chr_ids

  script:

"""
generate_multimappability_signal.csh ${fasta} ${params.read_size} ./
bash temp_map_reads.csh
bash temp_process_mapping.csh
"""
}

process music_preprocessing {
  tag "${file_id}"

  input:
    set file_id, file(bam) from bam_files
    file chr_ids from chr_ids.collect()

  output:
    set file_id, "preprocessed/*.tar" into preprocessed_bam_files

  script:

"""
mkdir preprocessed
samtools view *.bam | \
MUSIC -preprocess SAM stdin preprocessed/
mkdir preprocessed/sorted
MUSIC -sort_reads preprocessed/ preprocessed/sorted/
mkdir preprocessed/dedup
MUSIC -remove_duplicates ./preprocessed/sorted 2 preprocessed/dedup/
cd preprocessed
tar -c -f ${file_id}.tar *
"""
}

preprocessed_bam_files_control = Channel.create()
preprocessed_bam_files_chip = Channel.create()
preprocessed_bam_files.choice(
  preprocessed_bam_files_control,
  preprocessed_bam_files_chip ) { a -> a[0] =~ /.*control.*/ ? 0 : 1 }

process music_computation {
  tag "${file_id}"
  publishDir "results/peak_calling/${file_id}", mode: 'copy'

  input:
    set file_id, file(control) from preprocessed_bam_files_chip
    set file_id_control, file(chip) from preprocessed_bam_files_control.collect()
    file mapp from mappability.collect()

  output:
    file "*" into music_output_forward
    file "*.bed" into peaks_forward

  script:

"""
mkdir mappability control chip
mv ${mapp} mappability/
tar -xf ${control} -C control/
tar -xf ${chip} -C chip/

MUSIC -get_per_win_p_vals_vs_FC -chip chip/ -control control/ \
  -l_win_step ${params.step_l} \
  -l_win_min ${params.min_l} -l_win_max ${params.max_l}
MUSIC -get_multiscale_punctate_ERs \
  -chip chip/ -control control/ -mapp mappability/ \
  -l_mapp ${params.read_size} -l_frag ${params.frag_size} -q_val 1 -l_p 0
ls -l
"""
}
