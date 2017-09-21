# Add reads group
# Done locally due to problems with PSMN, using picard 2.12.1. 
# Obtained on 20 sept 2017: 
# git clone https://github.com/broadinstitute/picard.git
#    cd picard/
#    ./gradlew shadowJar
# copy build/libs/picard.jar to bin

##### Args #####
# $1: input BAM
# $2: id, ex: HKCVHBBXX.1
# $3: lb, ex: male
# $4: pu, ex: HKCVHBBXX.1
# $5: sm, ex: MRDR6
# $6: output BAM without RG in header

java -jar bin/picard.jar AddOrReplaceReadGroups I=$1 O=$6 RGID=$2 RGLB=$3 RGPL=illumina RGPU=$4 RGSM=$5 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true


