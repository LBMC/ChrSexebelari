#!/bin/sh

# generate training set
ll *.gz | sed 's/.gz//g' | awk '{system("gunzip -c "$9".gz > ~/data/JU28_59vs17_SNP/data/samples/"$9)}'
cd ~/data/JU28_59vs17_SNP/data/samples/
~/scripts/fastq_sampler/fastq_sampler.py -i 100000 -f MR_350_clean_1.fastq -g MR_350_clean_2.fastq
~/scripts/fastq_sampler/fastq_sampler.py -i 100000 -f MR_550_clean_1.fastq -g MR_550_clean_2.fastq
~/scripts/fastq_sampler/fastq_sampler.py -i 100000 -f NG-10944_JU2859_bis_lib169352_5217_1_1.fastq -g NG-10944_JU2859_bis_lib169352_5217_1_2.fastq
ll s_*.fastq | awk '{system("gzip < "$9" > "$9".gz")}'
cd ~/projects/JU28_59vs17_SNP/

# training set analysis

mkdir tests
cd tests
../nextflow ../src/SNP_calling.nf -c ../src/SNP_calling.config -profile docker --fasta "../data/fasta/DBG2OLC_output2.fasta" --fastq "../data/samples/*_{1,2}.fastq.gz" -resume -w ~/data/work_s/ --tumor "[\"s_NG-10944_JU2859_bis_lib169352_5217_1\"]" --normal "[\"s_MR_550_clean\", \"s_MR_350_clean\"]" --seq_number 800000
~/scripts/sms.sh "SNP done"

# real set analysis

./nextflow src/SNP_calling.nf -c src/SNP_calling.config -profile docker --fasta "data/fasta/DBG2OLC_output2.fasta" --fastq "data/fastq/*_{1,2}.fastq.gz" -resume -w ~/data/work/ --tumor "[\"NG-10944_JU2859_bis_lib169352_5217_1\"]" --normal "[\"MR_550_clean\", \"MR_350_clean\"]"
~/scripts/sms.sh "SNP done"
src/intersect_SNP.R \
  results/SNP/vcf_samtools/normal_sample_raw.csv \
  results/SNP/vcf_samtools/tumor_sample_filtered.csv \
  results/fasta/DBG2OLC_output2_filtered.fasta \
  data/list_of_enzymes.csv
~/scripts/sms.sh "SNP analysis done"

./nextflow src/SNP_calling.nf -c src/SNP_calling.config -profile docker --fasta "data/fasta/final_assembly.fasta" --fastq "data/fastq/*_{1,2}.fastq.gz" -resume -w ~/data/work/ --tumor "[\"NG-10944_JU2859_bis_lib169352_5217_1\"]" --normal "[\"MR_550_clean\", \"MR_350_clean\"]"
~/scripts/sms.sh "SNP done"
src/intersect_SNP.R \
  results/SNP/vcf_samtools/normal_sample_raw.csv \
  results/SNP/vcf_samtools/tumor_sample_filtered.csv \
  results/fasta/final_assembly_filtered.fasta \
  data/list_of_enzymes.csv
~/scripts/sms.sh "SNP analysis done"


# on the PSMN
find ~/data/ -name "MR_350_clean*"
find ~/data/ -name "MR_550_clean*"
find ~/data/ -name "10944_JU2859_bis_lib169352_5217_1*"

./nextflow src/SNP_calling.nf -c src/SNP_calling.config -profile sge --fasta "data/fasta/final_assembly.fasta" --fastq "/Xnfs/lbmcdb/Delattre_team/Request/Clean_data_for_assembly/*_{1,2}.fastq.gz" -resume -w /scratch/lmodolo/work/ --tumor "[\"NG-10944_JU2859_bis_lib169352_5217_1\"]" --normal "[\"MR_550_clean\", \"MR_350_clean\"]"

mkdir -p results/blastall/
makeblastdb -in data/fasta/DBG2OLC_output2.fasta -parse_seqids -dbtype nucl
blastn -query data/RNA5S_belari.fasta -db data/fasta/DBG2OLC_output2.fasta -out results/blastall/RNA5S_2.out
less results/blastall/RNA5S_2.out

makeblastdb -in data/fasta/DBG2OLC_output1.fasta -parse_seqids -dbtype nucl
blastn -query data/RNA5S_belari.fasta -db data/fasta/DBG2OLC_output1.fasta -out results/blastall/RNA5S_1.out
less results/blastall/RNA5S_1.out

makeblastdb -in data/fasta/nanoport_denovo.fasta -parse_seqids -dbtype nucl
blastn -query data/RNA5S_belari.fasta -db data/fasta/nanoport_denovo.fasta -out results/blastall/RNA5S_nanoport.out
less results/blastall/RNA5S_nanoport.out
