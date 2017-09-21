# Create sequence dictionnary
# Done locally due to problems with PSMN, using picard 2.12.1. picard.jar
# Obtained on 20 sept 2017: 
# git clone https://github.com/broadinstitute/picard.git
#    cd picard/
#    ./gradlew shadowJar
# copy build/libs/picard.jar to bin

##### Args #####
# $1: input ref in fasta
# $2: output dict 

java -jar bin/picard.jar CreateSequenceDictionary R=$1 O=$2 &&\

samtools faidx $1


