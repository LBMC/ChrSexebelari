export BUSCO_CONFIG_FILE=~/Programs/busco-master/config/config.ini
export PATH=~/Programs/augustus-3.3.1/bin/:$PATH
export PATH=~/Programs/augustus-3.3.1/scripts/:$PATH
export AUGUSTUS_CONFIG_PATH=~/Programs/augustus-3.3.1/config/

#SEQUENCE_file=~/Documents/stage_mbelari/results/hybrid_test/DBG2OLC/2018-06-01/final_assembly.fasta
SEQUENCE_file=~/Documents/stage_mbelari/results/annotation/BRAKER/soft_masked/augustus.hints.aa
OUTPUT_name=mbelari_prot
LINEAGE=~/Programs/busco-master/datasets/nematoda_odb9
MODE=proteins

python ~/Programs/busco-master/scripts/run_BUSCO.py -i ${SEQUENCE_file} -o ${OUTPUT_name} -l ${LINEAGE} -m ${MODE}
