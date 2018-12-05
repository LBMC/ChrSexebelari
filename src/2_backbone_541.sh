# compute backbone length

docker run -ti --volume "$(pwd)"/results:/root/results/ --volume "$(pwd)"/data:/root/data/ bioawk:1.0
bioawk -c fastx '{ print $name", "length($seq) }' \
  /root/data/fasta/final_assembly.fasta > /root/results/backbone_size.csv
exit

# get Backbone_541
docker run -ti --volume "$(pwd)"/results:/root/results/ --volume "$(pwd)"/data:/root/data/ bioawk:1.0
bioawk -c fastx '{ if($name = "Backbone_541") { print ">"$name; print $seq }}' \
  /root/data/fasta/final_assembly.fasta > /root/results/backbone_541.fasta
exit

# map read on Backbone_541
docker run -ti --volume "$(pwd)"/results:/root/results/ --volume "$(pwd)"/data:/root/data/ bowtie2:2.3.4.1
bowtie2-build --threads 8 /root/results/backbone_541.fasta /root/results/backbone_541.index

bowtie2 --very-sensitive -p 8 -x /root/results/backbone_541.index \
--rg-id MR_550 \
--rg PL:Illumina \
--rg SM:MR_550 \
-1 /root/results/fastq/trimming/MR_550_clean_trim_R1.fastq.gz \
-2 /root/results/fastq/trimming/MR_550_clean_trim_R2.fastq.gz | \
samblaster --addMateTags -M -i /dev/stdin | \
sambamba view -t 8 --valid -S -f bam -l 0 /dev/stdin \
-o /root/results/backbone_541_MR_550.bam

bowtie2 --very-sensitive -p 8 -x /root/results/backbone_541.index \
--rg-id MR_350 \
--rg PL:Illumina \
--rg SM:MR_350 \
-1 /root/results/fastq/trimming/MR_350_clean_trim_R1.fastq.gz \
-2 /root/results/fastq/trimming/MR_350_clean_trim_R2.fastq.gz | \
samblaster --addMateTags -M -i /dev/stdin | \
sambamba view -t 8 --valid -S -f bam -l 0 /dev/stdin \
-o /root/results/backbone_541_MR_350.bam

bowtie2 --very-sensitive -p 8 -x /root/results/backbone_541.index \
--rg-id MR_350 \
--rg PL:Illumina \
--rg SM:MR_350 \
-1 /root/results/fastq/trimming/NG-10944_JU2859_bis_lib169352_5217_1_trim_R1.fastq.gz \
-2 /root/results/fastq/trimming/NG-10944_JU2859_bis_lib169352_5217_1_trim_R2.fastq.gz | \
samblaster --addMateTags -M -i /dev/stdin | \
sambamba view -t 8 --valid -S -f bam -l 0 /dev/stdin \
-o /root/results/backbone_541_JU2859.bam
exit

