GENOME=~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta
BAM1=~/Documents/stage_mbelari/results/13genes/2018-06-12-MRDR5_vs_Hybrid_assembly.sorted.filtered.bam
BAM2=~/Documents/stage_mbelari/results/13genes/2018-06-12-MRDR6_vs_Hybrid_assembly.sorted.filtered.bam
OUTPUT=2018-07-24-MRDR_vs_Hybrid_assembly

#samtools mpileup -g -f ${GENOME} ${BAM1} ${BAM2} > ${OUTPUT}

bcftools view -v -c -g -O u -o ${OUTPUT}.var.bcf ${OUTPUT}.bcf 

#bcftools view ${OUTPUT}.var.bcf | vcfutils.pl varFilter - > ${OUTPUT}.var.final.vcf
