#$ -S /bin/bash
### change logs folder
#$ -o /home/cburny/logs
#$ -e /home/cburny/logs

module load nextflow/0.25.1

nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_sge --fastq_files "results/mapping/mapped/*.fastq.gz"
