#!/bin/bash

if [ $# -ne 1 ]
  then
    echo "usage: $0 <ligated_prefix>"
    exit -1
fi


BIN="sudo docker run --rm -v .:/data glimpse2-dev /bin/"
DATA=/data
PREFIX=$1

LST=ligate/$PREFIX.lst
ls -1v impute/imputed_*.bcf | sed 's/^/\/data\//g' > ${LST}

OUT=ligate/${PREFIX}_ligated.bcf
$BIN/GLIMPSE2_ligate \
	--input $DATA/${LST} \
	--output $DATA/$OUT
