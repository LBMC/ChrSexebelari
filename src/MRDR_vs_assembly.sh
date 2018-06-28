
assembly=~/Documents/stage_mbelari/results/redundans/2018-06-26-redundans_1/redundans/test1/scaffolds.reduced.fa
pe5_1=~/Documents/stage_mbelari/data/raw/2017_09_14_MRDR5_trim_Mbelari_mapped_sort_1.fastq
pe5_2=~/Documents/stage_mbelari/data/raw/2017_09_14_MRDR5_trim_Mbelari_mapped_sort_2.fastq
pe6_1=~/Documents/stage_mbelari/data/raw/2017_09_14_MRDR6_trim_Mbelari_mapped_sort_1.fastq
pe6_2=~/Documents/stage_mbelari/data/raw/2017_09_14_MRDR6_trim_Mbelari_mapped_sort_2.fastq
index_name=hybrid_reduced
out_name=2018-06-26-reduced_assembly_default

bowtie2-build ${assembly} ${index_name}

bowtie2 --very-sensitive -x ${index_name} -1 ${pe5_1} -2 ${pe5_2} -S ${out_name}_MRDR5.sam
bowtie2 --very-sensitive -x ${index_name} -1 ${pe6_1} -2 ${pe6_2} -S ${out_name}_MRDR6.sam

samtools view -Sb ${out_name}_MRDR5.sam > ${out_name}_MRDR5.bam
samtools view -Sb ${out_name}_MRDR6.sam > ${out_name}_MRDR6.bam

samtools sort ${out_name}_MRDR5.bam -o ${out_name}_MRDR5.sorted.bam
samtools sort ${out_name}_MRDR6.bam -o ${out_name}_MRDR6.sorted.bam

samtools index ${out_name}_MRDR5.sorted.bam
samtools index ${out_name}_MRDR6.sorted.bam


