##### Convert gff to bed : run locally
# $1: input gff
# $2: output bed

gff2bed < $1 > $2 &&\

bash src/date.sh $2
